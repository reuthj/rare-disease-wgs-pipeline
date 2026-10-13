#!/usr/bin/env python3
"""Acquire, prepare and validate the catalog without running analysis stages."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import urllib.request
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path, algorithm='sha256'):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def run(*args, stdout=None):
    result = subprocess.run([str(x) for x in args], stdout=stdout or subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode:
        raise RuntimeError(f'{args[0]} failed: {result.stderr.decode(errors="replace")[-3000:]}')
    return result.stdout.decode(errors='replace') if stdout is None else ''


def conda_executable():
    executable = os.environ.get('CONDA_EXE') or shutil.which('conda')
    if not executable:
        raise RuntimeError('Conda is required to run the predictor environment')
    return executable


def nonempty(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing/empty file: {path}')


def open_download(request, timeout):
    """Retry transient connection/server errors without relaxing integrity QC."""
    max_retries = 5
    for attempt in range(max_retries + 1):
        try:
            return urllib.request.urlopen(request, timeout=timeout)
        except HTTPError as exc:
            if exc.code not in (408, 429, 500, 502, 503, 504) or attempt == max_retries:
                raise
        except (URLError, ConnectionError, TimeoutError):
            if attempt == max_retries:
                raise
        delay = min(5 * (2 ** attempt), 60)
        print(f'Temporary download failure: {request.full_url}; retry {attempt + 1}/{max_retries} in {delay}s', flush=True)
        time.sleep(delay)


def text_open(path):
    return gzip.open(path, 'rt') if path.suffix in ('.gz', '.bgz') else path.open()


class RangeUnsupported(RuntimeError):
    """The server cannot provide safe segmented transfers."""


def parallel_download(source, partial, connections, timeout, chunk_bytes=8*1024*1024,
                      minimum_bytes=64*1024*1024):
    """Return transfer metadata, or None when sequential downloading is needed.

    A sidecar records source identity and hashes of completed ranges. File length
    alone never establishes completeness: the partial file may be sparse.
    """
    if connections <= 1 or not source.startswith(('https://', 'http://')):
        return None
    headers={'User-Agent':'rare-disease-resource-setup/1', 'Accept-Encoding':'identity'}
    try:
        with open_download(urllib.request.Request(source, headers=headers, method='HEAD'), timeout) as response:
            length=int(response.headers.get('Content-Length', '0'))
            etag=response.headers.get('ETag')
            if etag and etag.startswith('W/'):
                etag=None
            modified=response.headers.get('Last-Modified')
            if length < minimum_bytes or response.headers.get('Accept-Ranges','').lower() != 'bytes' or not (etag or modified):
                return None
            resolved=response.geturl()
    except HTTPError as exc:
        if exc.code in (405, 501):
            return None  # HEAD is explicitly unsupported; GET can still work.
        raise
    except OSError:
        # A transient connection failure does not establish that ranges are
        # unsupported. Preserve any partial transfer and its verified ledger.
        raise
    except ValueError:
        return None
    identity={'source':source, 'resolved_source':resolved, 'size':length,
              'etag':etag, 'last_modified':modified, 'chunk_bytes':chunk_bytes}
    ledger_path=Path(str(partial)+'.ranges.json')
    ledger={'identity':identity, 'completed':{}}
    try:
        previous=json.loads(ledger_path.read_text())
        if previous['identity']==identity and partial.stat().st_size==length:
            ledger=previous
        else:
            raise ValueError('Source or partial layout changed')
    except (OSError, ValueError, KeyError, TypeError):
        # A legacy sequential partial has no source ledger; preserve it for
        # inspection rather than treating its bytes as verified ranges.
        if partial.exists():
            os.replace(partial, Path(str(partial)+'.unverified-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')))
        ledger_path.unlink(missing_ok=True)
        with partial.open('wb') as stream:
            stream.truncate(length)
    count=(length+chunk_bytes-1)//chunk_bytes
    with partial.open('rb') as stream:
        for index, expected in list(ledger['completed'].items()):
            i=int(index)
            if not 0 <= i < count:
                del ledger['completed'][index]
                continue
            stream.seek(i*chunk_bytes)
            data=stream.read(min(chunk_bytes,length-i*chunk_bytes))
            if hashlib.sha256(data).hexdigest()!=expected:
                del ledger['completed'][index]
    lock=threading.Lock()
    stop=threading.Event()
    def save():
        temporary=Path(str(ledger_path)+'.tmp')
        temporary.write_text(json.dumps(ledger, sort_keys=True)+'\n')
        os.replace(temporary,ledger_path)
    save()
    pending=[i for i in range(count) if str(i) not in ledger['completed']]
    print(f'Downloading {partial.name}: {connections} connections; {count-len(pending)}/{count} verified ranges reused',flush=True)
    fd=os.open(partial,os.O_WRONLY)
    def worker():
        while not stop.is_set():
            with lock:
                if not pending:return
                i=pending.pop()
            start=i*chunk_bytes;end=min(length,start+chunk_bytes)-1
            request=urllib.request.Request(resolved,headers=dict(headers,Range=f'bytes={start}-{end}',**{'If-Range':etag or modified}))
            try:
                with open_download(request,timeout) as response:
                    if response.geturl()!=resolved:
                        raise ValueError('Segment source redirected after identity check')
                    if etag and response.headers.get('ETag')!=etag or not etag and response.headers.get('Last-Modified')!=modified:
                        raise ValueError('Source identity changed during segmented download')
                    if response.status==200:
                        raise RangeUnsupported('Server ignored byte ranges')
                    if response.status!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{length}':
                        raise ValueError('Incorrect HTTP Content-Range')
                    if response.headers.get('Content-Encoding','identity')!='identity':
                        raise ValueError('Unexpected HTTP content encoding')
                    data=response.read(end-start+2)
                    if len(data)!=end-start+1:
                        raise ValueError('Incomplete/oversized HTTP range')
                view=memoryview(data);offset=start
                while view:
                    written=os.pwrite(fd,view,offset)
                    if written<=0:raise OSError('Short write to partial download')
                    view=view[written:];offset+=written
                with lock:
                    ledger['completed'][str(i)]=hashlib.sha256(data).hexdigest()
                    save()
                    if len(ledger['completed']) % 32 == 0 or len(ledger['completed'])==count:
                        print(f'{partial.name}: {len(ledger["completed"])}/{count} ranges complete',flush=True)
            except Exception:
                stop.set()
                raise
    try:
        with ThreadPoolExecutor(max_workers=connections) as pool:
            futures=[pool.submit(worker) for _ in range(connections)]
            for future in futures:future.result()
        if len(ledger['completed'])!=count:
            raise ValueError('Missing ranges after download')
        os.fsync(fd)
    finally:
        os.close(fd)
    return dict(identity, download_date=now(), download_method='parallel_ranges', connections=connections)


def safe_archive_member(member, destination):
    parts=Path(member.name).parts
    if Path(member.name).is_absolute() or '..' in parts or member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
        raise ValueError(f'Unsafe archive member: {member.name}')
    return tarfile.data_filter(member,destination)


def validate(path, spec):
    """Lightweight format QC and pinned source identity; no reference comparisons."""
    nonempty(path)
    fmt = spec['format']
    for algorithm in ('sha256', 'md5'):
        if spec.get(algorithm) and digest(path, algorithm) != spec[algorithm]:
            raise ValueError(f'{algorithm} mismatch for preserved {spec["id"]}: {path}; '
                             'supply the matching snapshot or explicitly configure a new profile')
    if fmt == 'tar':
        with tarfile.open(path, 'r:gz') as archive:
            member=archive.next()
            if member is None:
                raise ValueError('Empty archive')
            safe_archive_member(member,path.parent)
        return
    if fmt.startswith('orphanet'):
        with path.open('rb') as stream:
            events=ET.iterparse(stream,events=('start','end'))
            _,root=next(events)
            if not root.attrib.get('date'):raise ValueError('Orphadata snapshot date missing')
            for event,node in events:
                if event!='end':continue
                if fmt=='orphanet6' and node.tag=='Gene' and node.findtext('Symbol'):return
                if fmt=='orphanet1' and node.tag=='ExternalReference' and node.findtext('Source','').strip().upper()=='MONDO' and node.findtext('Reference'):return
        raise ValueError(f'No usable mapping sample: {path}')
    if fmt == 'hpo':
        data = json.loads(path.read_text())
        if not data.get('graphs') or not data['graphs'][0].get('nodes'):
            raise ValueError('HPO graph/nodes absent')
        return
    if fmt == 'tbi':
        with gzip.open(path, 'rb') as stream:
            if stream.read(4) != b'TBI\x01':
                raise ValueError('Invalid Tabix index signature')
        return
    with text_open(path) as stream:
        lines = []
        for _ in range(150):
            line = stream.readline()
            if not line:
                break
            lines.append(line.rstrip('\n'))
    if not lines:
        raise ValueError('No readable records')
    body = [x for x in lines if x and not x.startswith('#')]
    if fmt == 'vcf':
        if not any(x.startswith('##fileformat=VCF') for x in lines):
            raise ValueError('VCF header missing')
        # bcftools reads the full header, which may exceed the preview.
        run('bcftools', 'view', '-h', path)
    elif fmt == 'gtf':
        if not body or len(body[0].split('\t')) != 9:
            raise ValueError('GTF requires nine columns')
        marker = spec.get('gtf_release_marker')
        if marker and marker.lower() not in '\n'.join(lines[:10]).lower():
            raise ValueError(f'GENCODE release metadata missing: {marker}')
    elif fmt == 'bed6':
        if not body or len(body[0].split('\t')) != 6:
            raise ValueError('Registry V4 requires six source columns')
        fields = body[0].split('\t')
        if not fields[4].startswith('EH38E') or int(fields[1]) >= int(fields[2]):
            raise ValueError('Unexpected Registry V4 ID/coordinate layout')
    elif fmt == 'bed4':
        if not body or len(body[0].split('\t')) != 4 or ':' not in body[0].split('\t')[3]:
            raise ValueError('Processed cCRE requires four columns and ID:CLASS')
    elif fmt == 'alphamissense':
        header = next((x.lstrip('#').split('\t') for x in lines if x.startswith('#CHROM\t')), [])
        if not {'CHROM','POS','REF','ALT','transcript_id','am_pathogenicity'}.issubset(header) or not body:
            raise ValueError('AlphaMissense required columns/records missing')
    elif fmt in ('gencc', 'constraint', 'hpoa'):
        header_line = next((x for x in lines if x.startswith('#database_id\t')), None) if fmt == 'hpoa' else None
        header = (header_line.lstrip('#') if header_line else body[0] if body else '').split('\t')
        requirements = {'gencc': [('gene_symbol','submitted_as_hgnc_symbol'), ('classification_title','submitted_as_classification_name'), ('disease_title','submitted_as_disease_name'), ('submitter_title','submitted_as_submitter_name')],
                        'constraint': [('gene','gene_symbol'), ('transcript',), ('mis.oe','mis_oe'), ('mis.z_score','mis_z'), ('lof.oe_ci.upper','lof_oe_ci_upper')],
                        'hpoa': [('database_id',), ('hpo_id',)]}[fmt]
        if any(not set(group).intersection(header) for group in requirements):
            raise ValueError(f'{fmt} required columns missing: {header}')
        if not body:
            raise ValueError('No data rows')
    elif fmt == 'acmg':
        rows = list(csv.reader(lines))
        if spec.get('report_marker', 'c3.3') not in lines[0] or len(rows)<3 or not {'Gene Symbol','MONDO ID','MOI','Reportable as SF'}.issubset(rows[1]):
            raise ValueError('Expected ACMG SF v3.3 CSV, not HTML or a different release')
    elif fmt == 'perl' and 'package NMD' not in '\n'.join(lines):
        raise ValueError('NMD Perl package missing')


class Bootstrap:
    def __init__(self, root, catalog, args):
        self.root, self.catalog, self.args = root, catalog, args
        self.state_path = root / 'data/resources/resource_manifest.json'
        self.state = json.loads(self.state_path.read_text()) if self.state_path.exists() else {'schema_version': 1, 'resources': {}}
        self.specs = {r['id']: r for r in catalog['resources']}
        self.complete = set()
        for r in catalog['resources']:
            if r['id'] in self.state['resources']:
                continue
            status='NOT_CHECKED'
            if r.get('processed'):
                processed=self.path(r['processed']);raw=self.path(r['raw'])
                status='UNVALIDATED' if processed.exists() else 'DOWNLOADED_NOT_PROCESSED' if raw.exists() else 'PARTIAL' if Path(str(raw)+'.partial').exists() else 'MISSING'
            self.state['resources'][r['id']]={k:r.get(k) for k in ('id','version','sources','raw','processed','required')}
            self.state['resources'][r['id']]['status']=status

    def path(self, relative):
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root):
            raise ValueError(f'Destination outside project: {relative}')
        return path

    def predictor_python_command(self):
        """Use the predictor beside this main environment by exact prefix."""
        sibling = Path(sys.prefix).resolve().parent / self.args.predictor_env
        if sibling.is_dir():
            return [conda_executable(), 'run', '--no-capture-output', '--prefix', str(sibling), 'python']
        if Path(sys.prefix).resolve().parent == (self.root / '.conda/envs').resolve():
            raise RuntimeError(f'Project-local predictor environment is missing: {sibling}; create it from envs/spliceai-env.yaml')
        return [conda_executable(), 'run', '--no-capture-output', '-n', self.args.predictor_env, 'python']

    def save(self):
        self.state['checked_at'] = now()
        self.state['profile'] = self.catalog.get('profile')
        self.state['catalog_sha256'] = hashlib.sha256(json.dumps(self.catalog, sort_keys=True).encode()).hexdigest()
        self.state['limitations'] = self.catalog.get('limitations', [])
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.state_path.with_suffix('.json.partial')
        tmp.write_text(json.dumps(self.state, indent=2) + '\n')
        os.replace(tmp, self.state_path)

    def recipe_key(self, r):
        deps = {d: self.state['resources'].get(d, {}).get('artifacts', {}) for d in r.get('dependencies', [])}
        return hashlib.sha256(json.dumps([r, deps], sort_keys=True).encode()).hexdigest()

    def file_record(self, path):
        stat = path.stat()
        return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'sha256': digest(path)}

    def outputs(self, r):
        p = self.path(r['processed'])
        return [p] + [Path(str(p) + ext) for ext in r.get('indexes', [])]

    def cache_inventory(self, r):
        base = self.path(r['processed'])
        entries = [(str(p.relative_to(base)), p.stat().st_size, p.stat().st_mtime_ns)
                   for p in sorted(base.rglob('*')) if p.is_file()]
        return {'files': len(entries), 'total_bytes': sum(x[1] for x in entries),
                'inventory_sha256': hashlib.sha256(json.dumps(entries).encode()).hexdigest()}

    def prepared_valid(self, r):
        p = self.path(r['processed'])
        if r['kind'] == 'vep_cache':
            info = p / 'info.txt'
            nonempty(info)
            text = info.read_text()
            if r.get('cache_gencode', 'GENCODE50') not in text or r.get('assembly', 'GRCh38') not in text:
                raise ValueError('VEP cache build/GENCODE mismatch')
            for chrom in [str(i) for i in range(1,23)] + ['X','Y']:
                chunks = list((p / chrom).glob('*.gz'))
                if not chunks:
                    raise ValueError(f'Incomplete VEP cache chromosome {chrom}')
                nonempty(chunks[0])
                with gzip.open(chunks[0], 'rb') as stream:stream.read(100)
            return
        # Source checksums were enforced during acquisition. Routine output QC
        # samples structure; recorded size/mtime handles inexpensive READY reuse.
        processed = dict(r,sha256=None,md5=None)
        if r['kind'] == 'fasta':
            nonempty(p)
            with p.open() as stream:
                if not stream.readline().startswith('>1'):
                    raise ValueError('Expected no-chr GRCh38 FASTA')
            nonempty(Path(str(p)+'.fai'))
            names = {x.split('\t')[0] for x in Path(str(p)+'.fai').read_text().splitlines()}
            if not {'1','X','Y','M'}.issubset(names) and not self.args.fixture:
                raise ValueError('FASTA primary chromosome convention mismatch')
            run('samtools','faidx',p,'1:1-10')
            index = Path(str(p)+'.fxi'); nonempty(index)
            with sqlite3.connect(f'file:{index}?mode=ro', uri=True) as db:
                if not db.execute('SELECT name FROM sqlite_master LIMIT 1').fetchone():
                    raise ValueError('Empty/unreadable pyfastx index')
            return
        if r['kind'] == 'ccre':
            processed.update(format='bed4', sha256=None, md5=None)
        elif r['kind'] == 'clinvar':
            processed.update(sha256=None, md5=None)
        validate(p, processed)
        for index in self.outputs(r)[1:]:
            nonempty(index)
        if r.get('indexes') == ['.tbi']:
            contigs = run('tabix','-l',p).splitlines()
            if not contigs:
                raise ValueError('Empty Tabix contig list')
            run('tabix',p,contigs[0]+':1-100000')
            if r['kind'] == 'clinvar' and any(c not in [str(i) for i in range(1,23)]+['X','Y','M'] for c in contigs):
                raise ValueError('Normalized ClinVar has non-primary/incorrect contigs')

    def reusable(self, r, previous):
        if previous.get('recipe_key') != self.recipe_key(r) or previous.get('status') != 'READY':
            return False
        try:
            for rel, identity in previous.get('artifacts', {}).items():
                p = self.path(rel); stat = p.stat()
                if stat.st_size != identity['size'] or stat.st_mtime_ns != identity['mtime_ns']:
                    return False
            if not previous.get('artifacts'):
                return False
            self.prepared_valid(r)
            if r['kind'] == 'vep_cache' and previous.get('cache_inventory') != self.cache_inventory(r):
                return False
            return True
        except (OSError, ValueError, RuntimeError, sqlite3.Error):
            return False

    def acquire_bundled_snapshot(self, r, dest, spec):
        """Import a pinned bundled export without extracting arbitrary members."""
        folder = self.root / 'config/resource_snapshots'
        manifest = json.loads((folder / 'manifest.json').read_text())
        entry = next(item for item in manifest['snapshots'] if item['id'] == r['id'])
        if entry['sha256'] != r.get('sha256') or entry['file'] != dest.name:
            raise ValueError(f'Bundled snapshot does not match catalog: {r["id"]}')
        archive = folder / manifest['archive']
        if archive.parent.resolve() != folder.resolve() or digest(archive) != manifest['archive_sha256']:
            raise ValueError('Bundled snapshot archive path/checksum mismatch')
        with tempfile.TemporaryDirectory(prefix='.bundled-snapshot-', dir=dest.parent) as temp:
            staged = Path(temp) / dest.name
            with tarfile.open(archive, 'r:gz') as bundle:
                members = [member for member in bundle.getmembers() if member.name == entry['file']]
                if len(members) != 1 or not members[0].isfile() or members[0].size != entry['bytes']:
                    raise ValueError('Missing, duplicate or invalid bundled snapshot member')
                with bundle.extractfile(members[0]) as stream, staged.open('wb') as target:
                    shutil.copyfileobj(stream, target)
            validate(staged, spec)
            os.replace(staged, dest)
        return dest, {'acquisition': 'bundled_snapshot', 'source': str(archive.relative_to(self.root)),
                      'archive_sha256': manifest['archive_sha256'], 'snapshot_sha256': entry['sha256'],
                      'download_date': None, 'acquired_at': now()}

    def acquire(self, r, url=None, destination=None, spec=None):
        dest = destination or self.path(r['raw'])
        spec = spec or r
        previous = self.state['resources'].get(r['id'], {})
        recorded = previous.get('artifacts', {}).get(str(dest.relative_to(self.root)))
        if not recorded and previous.get('archive_identity', {}).get('path') == str(dest.relative_to(self.root)):
            recorded = previous['archive_identity']
        current_spec=hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest()
        artifact_spec=previous.get('artifact_spec_sha256',previous.get('spec_sha256'))
        changed_unpinned=recorded and artifact_spec!=current_spec and not spec.get('sha256') and not spec.get('md5')
        if recorded and artifact_spec == current_spec:
            spec = dict(spec, sha256=spec.get('sha256') or recorded['sha256'])
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            if changed_unpinned:
                raise ValueError('Explicit new unpinned specification requires fresh acquisition')
            validate(dest, spec)
            return dest, {'acquisition': 'existing_valid_file', 'download_date': None, 'source': None}
        except (OSError, ValueError, EOFError, tarfile.TarError, RuntimeError):
            pass
        if (url is None and self.args.snapshot_dir is None
                and self.catalog.get('profile') == 'preserved-comparison'
                and r['id'] in ('gencc', 'acmg_sf')
                and (self.root / 'config/resource_snapshots').exists()):
            return self.acquire_bundled_snapshot(r, dest, spec)
        sources = [url] if url else r['sources']
        if self.args.snapshot_dir:
            snapshot = self.args.snapshot_dir / dest.name
            if snapshot.exists():
                sources = [snapshot.resolve().as_uri()] + sources
        partial = Path(str(dest) + '.partial')
        errors = []
        for source in sources:
            try:
                try:
                    metadata=parallel_download(source,partial,self.args.download_connections,self.args.timeout)
                except RangeUnsupported as exc:
                    print(f'{dest.name}: {exc}; falling back to one connection',flush=True)
                    metadata=None
                if metadata is None:
                    request = urllib.request.Request(source, headers={'User-Agent':'rare-disease-resource-setup/1', 'Accept-Encoding':'identity'})
                    with open_download(request, self.args.timeout) as response, partial.open('wb') as stream:
                        Path(str(partial)+'.ranges.json').unlink(missing_ok=True)
                        shutil.copyfileobj(response, stream, length=4*1024*1024)
                        expected = response.headers.get('Content-Length')
                        metadata = {'source': source, 'resolved_source': response.geturl(), 'download_date': now(),
                                    'etag': response.headers.get('ETag'), 'last_modified': response.headers.get('Last-Modified'),
                                    'download_method':'sequential'}
                        if source.startswith('file:'):
                            metadata.update(acquisition='supplied_snapshot', download_date=None, acquired_at=now())
                    if expected and partial.stat().st_size != int(expected):
                        raise ValueError('Incomplete Content-Length')
                # Keep original extension for compression-aware validation.
                probe = dict(spec)
                # Validation relies on extension for text resources: temporary staging symlink is unnecessary.
                staging = dest.parent / ('.validate-' + dest.name)
                os.replace(partial, staging)
                try:
                    validate(staging, probe)
                    os.replace(staging, dest)
                    Path(str(partial)+'.ranges.json').unlink(missing_ok=True)
                finally:
                    if staging.exists(): staging.unlink()
                return dest, metadata
            except (OSError, ValueError, EOFError, tarfile.TarError, RuntimeError) as exc:
                errors.append(f'{source}: {exc}')
        raise RuntimeError('Download/integrity failed: ' + '; '.join(errors))

    def prepare(self, r):
        # A structurally valid installed cache can be adopted without reacquiring
        # a multi-GB archive. Its release/build and complete file inventory are recorded.
        if r['kind'] == 'vep_cache' and not self.path(r['raw']).exists():
            try:
                self.prepared_valid(r)
                previous = self.state['resources'].get(r['id'], {})
                if previous.get('artifact_spec_sha256', previous.get('spec_sha256')) not in (None, hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest()):
                    raise ValueError('Changed VEP specification requires acquisition')
                if previous.get('cache_inventory') and previous['cache_inventory'] != self.cache_inventory(r):
                    raise ValueError('Previously prepared VEP cache tree changed; reacquire archive')
                if previous.get('archive_identity'):
                    return previous['acquisition']
                return {'acquisition':'adopted_installed_cache', 'download_date':None,
                        'source':None, 'identity_limit':'Historical archive checksum unavailable; release/build and installed tree inventoried'}
            except (OSError, ValueError):
                pass
        raw, acquisition = self.acquire(r)
        dest = self.path(r['processed']); dest.parent.mkdir(parents=True, exist_ok=True)
        if r['kind']=='copy' and raw!=dest:
            with tempfile.TemporaryDirectory(prefix='.resource-copy-',dir=dest.parent) as temp:
                staged=Path(temp)/dest.name
                shutil.copyfile(raw,staged)
                validate(staged,r)
                os.replace(staged,dest)
        if r.get('md5_source'):
            md5path, _ = self.acquire(r, r['sources'][0]+'.md5', Path(str(raw)+'.md5'),
                                      {'id':r['id'], 'format':'text'})
            expected = md5path.read_text().split()[0]
            if len(expected) != 32 or digest(raw,'md5') != expected:
                raise ValueError('ClinVar upstream MD5 mismatch')
        if r.get('raw_index'):
            self.acquire(r, r['sources'][0]+'.tbi', Path(str(raw)+'.tbi'), {'id':r['id'],'format':'tbi'})
        if r.get('supplied_index'):
            self.acquire(r, r['sources'][0]+'.tbi', Path(str(dest)+'.tbi'), {'id':r['id'],'format':'tbi'})
        elif r['kind'] == 'copy' and r.get('indexes'):
            run('tabix','-f','-s','1','-b','2','-e','2','-S','1',dest)
        elif r['kind'] in ('fasta','ccre','clinvar','vep_cache'):
            with tempfile.TemporaryDirectory(prefix='.resource-build-',dir=dest.parent) as temp:
                tmp = Path(temp)
                if r['kind'] == 'fasta':
                    staged = tmp / dest.name
                    with gzip.open(raw,'rt') as source, staged.open('w') as target:
                        for line in source:
                            target.write('>'+line[4:] if line.startswith('>chr') else line)
                    run('samtools','faidx',staged)
                    command = [self.args.pyfastx_python] if self.args.pyfastx_python else self.predictor_python_command()
                    run(*command,'-c','import pyfastx,sys; pyfastx.Fasta(sys.argv[1])',staged)
                    for ext in ('.fai','.fxi'): os.replace(Path(str(staged)+ext),Path(str(dest)+ext))
                    os.replace(staged,dest)
                elif r['kind'] == 'ccre':
                    bed = tmp / 'ccre.bed'
                    with raw.open() as source, bed.open('w') as target:
                        for line in source:
                            if not line.strip() or line.startswith('#'): continue
                            f = line.rstrip('\n').split('\t')
                            if len(f)!=6 or not f[4].startswith('EH38E') or int(f[1])>=int(f[2]):
                                raise ValueError('Unexpected cCRE source layout')
                            target.write('\t'.join(f[:3]+[f[4]+':'+f[5]])+'\n')
                    run('sort','-k1,1V','-k2,2n',bed,'-o',bed)
                    staged = tmp / dest.name
                    with staged.open('wb') as target: run('bgzip','-c',bed,stdout=target)
                    run('tabix','-p','bed',staged)
                    os.replace(Path(str(staged)+'.tbi'),Path(str(dest)+'.tbi')); os.replace(staged,dest)
                elif r['kind'] == 'clinvar':
                    primary, renamed, staged = tmp/'primary.vcf.gz',tmp/'renamed.vcf.gz',tmp/dest.name
                    run('bcftools','view','-r',','.join([str(i) for i in range(1,23)]+['X','Y','MT']),'-Oz','-o',primary,raw)
                    mapping=tmp/'rename.txt';mapping.write_text('MT\tM\n')
                    run('bcftools','annotate','--rename-chrs',mapping,'-Oz','-o',renamed,primary)
                    reference = self.path(self.specs['reference']['processed'])
                    run('bcftools','norm','-f',reference,'-c','e','-m','-any','-Oz','-o',staged,renamed)
                    run('tabix','-p','vcf',staged)
                    os.replace(Path(str(staged)+'.tbi'),Path(str(dest)+'.tbi'));os.replace(staged,dest)
                else:
                    # Reject unsafe members during extraction, without a separate archive scan.
                    with tarfile.open(raw,'r:gz') as archive: archive.extractall(tmp,filter=safe_archive_member)
                    staged=tmp/dest.parent.name/dest.name
                    if not (staged/'info.txt').is_file(): raise ValueError('VEP archive layout mismatch')
                    if dest.exists():
                        quarantine=dest.with_name(dest.name+'.incomplete-'+datetime.now().strftime('%Y%m%d%H%M%S'))
                        os.replace(dest,quarantine)
                    os.replace(staged,dest)
        self.prepared_valid(r)
        return acquisition

    def check_software_assets(self, r):
        code = ('import importlib.metadata as m, pathlib, spliceai; '
                'assert m.version("spliceai")=="1.3.1"; '
                'p=pathlib.Path(spliceai.__file__).parent; '
                'models=list(p.rglob("*.h5")); assert len(models)>=5; '
                'print("models",len(models))')
        return {'probe':run(*self.predictor_python_command(),'-c',code).strip()}

    def execute(self, id):
        if id in self.complete: return
        r=self.specs[id]
        for dep in r.get('dependencies',[]): self.execute(dep)
        previous=self.state['resources'].get(id,{})
        record={k:r.get(k) for k in ('id','version','sources','raw','processed','indexes','consumers','required','limitation')}
        record['checked_at']=now()
        record['spec_sha256']=hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest()
        try:
            if any(self.state['resources'][d]['status']!='READY' for d in r.get('dependencies',[])):
                raise ValueError('Required preparation dependency is not READY')
            if r['kind'].startswith('remote_'):
                record.update(status='DEFERRED_TO_ANALYSIS_STAGE', disposition='Remote access occurs only in its consuming analysis stage; no setup probe performed')
            elif r['kind']=='software_assets':
                record.update(self.check_software_assets(r))
            else:
                if self.reusable(r,previous):
                    record.update({k:previous[k] for k in ('artifacts','recipe_key','acquisition','cache_inventory','artifact_spec_sha256','archive_identity','archive_cleanup') if k in previous})
                    record['action']='REUSED'
                elif self.args.check_only:
                    raise ValueError('No valid prepared manifest entry; run setup to validate/adopt or rebuild resource')
                else:
                    record['status']='PREPARING'
                    self.state['resources'][id]=dict(previous, status='PREPARING')
                    self.save()
                    acquisition=self.prepare(r)
                    record.update(acquisition=acquisition,recipe_key=self.recipe_key(r),action='PREPARED')
                    paths=self.outputs(r) if r['kind']!='vep_cache' else [self.path(r['processed'])/'info.txt']
                    if self.path(r['raw']).is_file() and self.path(r['raw']) not in paths:paths.append(self.path(r['raw']))
                    for ext in ('.tbi','.md5') if r.get('raw_index') else ():
                        paths.append(Path(str(self.path(r['raw']))+ext))
                    record['artifacts']={str(p.relative_to(self.root)):self.file_record(p) for p in paths}
                    record['artifact_spec_sha256']=record['spec_sha256']
                    if r['kind']=='vep_cache':record['cache_inventory']=self.cache_inventory(r)
                    if r['kind']=='vep_cache' and previous.get('archive_identity') and not self.path(r['raw']).exists():
                        record['archive_identity']=previous['archive_identity']
            if not r['kind'].startswith('remote_'):record['status']='READY'
        except Exception as exc:
            record.update(status='INVALID_OR_INCOMPLETE',step='readiness' if self.args.check_only else 'acquisition/preparation/readiness',reason=str(exc),next_step='Correct source/artifact or prerequisites and rerun the same command; see docs/sources_setup.md')
            # Preserve prior identities even after failure, so a subsequent retry
            # cannot silently accept a changed live export under the same recipe.
            for key in ('artifacts','acquisition','recipe_key','cache_inventory','artifact_spec_sha256','archive_identity','archive_cleanup'):
                if key in previous: record[key]=previous[key]
            if previous.get('artifacts') and 'artifact_spec_sha256' not in record:
                record['artifact_spec_sha256']=previous.get('spec_sha256')
        self.state['resources'][id]=record;self.complete.add(id);self.save()
        print(id+': '+record['status']+' '+record.get('action',record.get('reason','')),flush=True)

    def cleanup_vep_archive(self, id):
        """Retain download identity durably before removing a validated cache archive."""
        r=self.specs[id];record=self.state['resources'][id]
        if r['kind']!='vep_cache' or record['status']!='READY':return
        raw=self.path(r['raw']);relative=str(raw.relative_to(self.root))
        identity=record.get('artifacts',{}).get(relative)
        if identity:
            record['archive_identity']=dict(identity,path=relative)
            del record['artifacts'][relative]
            record['archive_cleanup']='PENDING'
            self.save()
        if not record.get('archive_identity'):return
        # No acquisition or extraction here: execute() has already validated the
        # installed tree. A crash before/after unlink is safe on the next run.
        raw.unlink(missing_ok=True)
        record['archive_cleanup']='REMOVED_AFTER_CACHE_VALIDATION'
        self.save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--catalog',type=Path)
    parser.add_argument('--snapshot-dir',type=Path,help='Optional historical raw artifacts, matched by basename and checksum')
    parser.add_argument('--only',nargs='+',help='Subset for maintenance; never reports overall READY')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--timeout',type=int,default=60)
    parser.add_argument('--download-connections',type=int,default=4,help='Parallel connections for large HTTP downloads (1 disables; default: 4)')
    parser.add_argument('--predictor-env',default='spliceai-env')
    parser.add_argument('--pyfastx-python',help='Python interpreter with pyfastx, otherwise use predictor env')
    parser.add_argument('--fixture',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();root=args.project_root.resolve();root.mkdir(parents=True,exist_ok=True)
    if not 1 <= args.download_connections <= 16:parser.error('--download-connections must be between 1 and 16')
    catalog_path=args.catalog or root/'config/resource_catalog.json'
    catalog=json.loads(catalog_path.read_text())
    if args.snapshot_dir is None and catalog.get('profile')=='preserved-comparison' and (root/'data/comparison_snapshots').is_dir():
        args.snapshot_dir=root/'data/comparison_snapshots'
    if len({r['id'] for r in catalog['resources']})!=len(catalog['resources']):parser.error('Duplicate resource IDs')
    selected=args.only or [r['id'] for r in catalog['resources']]
    if not set(selected).issubset({r['id'] for r in catalog['resources']}):parser.error('Unknown --only resource')
    # Paths with scientific consumers must agree with the current workflow configuration.
    if (root/'config/config.yaml').exists():
        import yaml
        config=yaml.safe_load((root/'config/config.yaml').read_text())
        for r in catalog['resources']:
            if r.get('config_key'):
                value=config
                for key in r['config_key'].split('.'):value=value[key]
                expected=r['sources'][0] if r['kind'].startswith('remote_') else r['processed']
                if value!=expected:parser.error(f'Path/source mismatch for {r["id"]}: config {value!r}, catalog {expected!r}; update both consistently')
        for r in catalog['resources']:
            if r['id']=='clinvar':
                for key, expected in [('clinvar_vcf',r['raw']),('clinvar_index',r['raw']+'.tbi'),('clinvar_md5',r['raw']+'.md5')]:
                    if config['stage5']['resources'][key]!=expected:parser.error(f'ClinVar source path mismatch: stage5.resources.{key}')
            if r['id']=='vep_cache':
                vep=config['stage1']['vep']
                if str(vep['cache_version'])!=str(r['cache_version']) or vep['assembly']!=r['assembly'] or r['processed']!=vep['cache_dir']+'/homo_sapiens/'+str(r['cache_version'])+'_'+r['assembly']:
                    parser.error('VEP cache catalog/config release, assembly or directory mismatch')
                if str(vep['version']).split('.')[0]!=str(r['cache_version']):
                    parser.error('Configured VEP software and cache major releases differ')
                gencode=next((x for x in catalog['resources'] if x['id']=='gencode'),None)
                if gencode and r['cache_gencode'].replace(' ','')!='GENCODE'+str(gencode['version']):
                    parser.error('VEP cache and GENCODE catalog releases differ')
    if not args.fixture:
        from install_vep import install_vep
        try:
            install_vep(config['stage1']['vep']['version'] if (root/'config/config.yaml').exists() else '116.1',
                        check_only=args.check_only)
        except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
            parser.exit(1, f'VEP setup failed: {exc}\n')
    lock=root/'data/resources/.setup.lock';lock.parent.mkdir(parents=True,exist_ok=True)
    with lock.open('a') as stream:
        try:fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:parser.error('Resource setup already running for this project')
        bootstrap=Bootstrap(root,catalog,args)
        bootstrap.state['status']='CHECKING' if args.check_only else 'PREPARING'
        bootstrap.save()
        for id in selected:bootstrap.execute(id)
        if not args.check_only:
            for id in sorted(bootstrap.complete):
                try:bootstrap.cleanup_vep_archive(id)
                except OSError as exc:
                    bootstrap.state['resources'][id].update(status='INVALID_OR_INCOMPLETE',step='archive cleanup',reason=str(exc),next_step='Check archive permissions and rerun setup')
                    bootstrap.save()
        failed=[id for id in bootstrap.complete if bootstrap.specs[id].get('required',True) and not bootstrap.specs[id]['kind'].startswith('remote_') and bootstrap.state['resources'][id]['status']!='READY']
        status='INCOMPLETE' if failed else 'SUBSET_READY' if args.only else 'READY'
        bootstrap.state['status']=status;bootstrap.save();print(status,flush=True)
        return 1 if failed else 0


if __name__=='__main__':
    sys.exit(main())
