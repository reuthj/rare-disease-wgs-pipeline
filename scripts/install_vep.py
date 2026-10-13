#!/usr/bin/env python3
"""Install the workflow's pinned VEP software into the active main Conda environment."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

VERSION = '116.1'
# Official source revisions from the functioning VEP 116.1 installation.
SOURCES = {
    'ensembl-vep': '97c24b205f4ba8b44ccc75b839512a1d59104019',
    'ensembl': '0d852313c420a1d5128b6774635125ebbac03350',
    'ensembl-variation': 'db7394a02523a6e519cb0d89f84b2889bac12737',
    'ensembl-funcgen': '90049ea7ee4d8ae3a6d298dca46d6c6ab20538c4',
    'ensembl-compara': '8057d0ebb016bc4331dc6e071916f3fd594f5abd',
    'ensembl-io': '6afb5dc27a5ae6881d75959153fbe6e9a4a7e788',
    'Bio-DB-HTS': '5daa12aa3dd7651f664830d9e1e10cea1c469b2f',
}
HTS_CHECK = ['-MBio::DB::HTS', '-MBio::DB::HTS::Faidx',
             '-MBio::DB::HTS::Tabix', '-e', 'exit($Bio::DB::HTS::VERSION eq "3.01" ? 0 : 1)']


def software_environment(prefix):
    env = os.environ.copy()
    for key in ('PERL5LIB', 'PERL_LOCAL_LIB_ROOT', 'PERL_MB_OPT', 'PERL_MM_OPT'):
        env.pop(key, None)
    env['PATH'] = str(prefix / 'bin') + os.pathsep + env.get('PATH', '')
    # Perl's build tools expect the compiler driver, rather than raw ld.
    env['CC'] = str(prefix / 'bin/x86_64-conda-linux-gnu-cc')
    env['LD'] = env['CC']
    return env


def ready(prefix, env):
    try:
        help_text = subprocess.check_output([prefix / 'bin/vep', '--help'],
                                            env=env, stderr=subprocess.STDOUT, text=True)
        if not any(line.strip() == f'ensembl-vep          : {VERSION}'
                   for line in help_text.splitlines()):
            return False
        subprocess.run([prefix / 'bin/perl', *HTS_CHECK], env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (OSError, subprocess.CalledProcessError):
        return False


def install_vep(expected_version=VERSION, check_only=False):
    if str(expected_version) != VERSION:
        raise RuntimeError(f'This installer provides VEP {VERSION}; update software and cache versions together')
    prefix = Path(sys.prefix).resolve()
    if not (prefix / 'conda-meta').is_dir():
        raise RuntimeError('Activate the main project Conda environment before installing VEP')
    env = software_environment(prefix)
    if ready(prefix, env):
        print(f'VEP {VERSION}: READY (reused)', flush=True)
        return
    if check_only:
        raise RuntimeError(f'VEP {VERSION} is incomplete; rerun setup without --check-only')
    if not Path(env['CC']).is_file():
        raise RuntimeError('Main environment compiler is missing; recreate it from envs/environment.yaml')
    share = prefix / 'share'
    share.mkdir(exist_ok=True)
    with (share / '.vep-install.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if ready(prefix, env):
            return
        print(f'Installing VEP {VERSION} into the main environment', flush=True)
        with tempfile.TemporaryDirectory(prefix='vep-install-', dir=share) as temporary:
            staging = Path(temporary)
            trees, records = {}, {}
            for repo, revision in SOURCES.items():
                url = f'https://codeload.github.com/Ensembl/{repo}/tar.gz/{revision}'
                archive = staging / f'{repo}.tar.gz'
                print(f'  {repo}', flush=True)
                with urllib.request.urlopen(url, timeout=120) as response, archive.open('wb') as output:
                    shutil.copyfileobj(response, output)
                with archive.open('rb') as stream:
                    checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
                destination = staging / repo
                destination.mkdir()
                with tarfile.open(archive) as source:
                    source.extractall(destination, filter='data')
                trees[repo] = next(destination.iterdir())
                records[repo] = {'revision': revision, 'source': url, 'sha256': checksum}
            for command in ([prefix / 'bin/perl', 'Build.PL', '--htslib', prefix],
                            [prefix / 'bin/perl', 'Build'],
                            [prefix / 'bin/perl', 'Build', 'install']):
                subprocess.run(command, cwd=trees['Bio-DB-HTS'], env=env, check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            vep = trees['ensembl-vep']
            api = vep / 'Bio/EnsEMBL'
            api.mkdir(parents=True)
            for repo in SOURCES:
                if repo in ('ensembl-vep', 'Bio-DB-HTS'):
                    continue
                shutil.copytree(trees[repo] / 'modules/Bio/EnsEMBL', api, dirs_exist_ok=True)
                metadata = vep / '.version'
                metadata.mkdir(exist_ok=True)
                (metadata / repo).write_text(f'release 116\nsub {SOURCES[repo]}\n')
            (vep / 'software_sources.json').write_text(json.dumps(records, indent=2) + '\n')
            target = share / f'ensembl-vep-{VERSION}'
            if target.exists():
                shutil.rmtree(target)
            shutil.move(str(vep), target)
            launcher = prefix / 'bin/vep'
            launcher.write_text(
                '#!/bin/sh\n'
                'base=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)\n'
                'unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT\n'
                f'exec "$base/bin/perl" "$base/share/ensembl-vep-{VERSION}/vep" "$@"\n')
            launcher.chmod(0o755)
        if not ready(prefix, env):
            raise RuntimeError('VEP installation failed its startup check; recreate the main environment and rerun setup')
    print(f'VEP {VERSION}: READY', flush=True)


if __name__ == '__main__':
    try:
        install_vep()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        if isinstance(error, subprocess.CalledProcessError) and error.stdout:
            print(error.stdout.decode(errors='replace')[-3000:], file=sys.stderr)
        print(f'VEP setup failed: {error}', file=sys.stderr)
        sys.exit(1)
