# Stage 6 Scoring Specification

## Purpose

Stage 6 prioritizes retained rare-disease WGS candidates using a
transparent **0--100 final score**. The score combines five independent
evidence components:

  Component                   Weight in Final Score
  ------------------------- -----------------------
  Variant pathogenicity                         30%
  Mode of inheritance                           20%
  Gene / disease evidence                       25%
  Case relevance                                15%
  Rarity                                        10%

The final score is:

**Final Score = 0.30 × Variant Pathogenicity + 0.20 × Inheritance + 0.25
× Gene/Disease Evidence + 0.15 × Case Relevance + 0.10 × Rarity**

Each component is scored on a **0--100 scale** before weighting.

For evidence measures that use continuous numerical predictors,
predefined anchor points are used and values between anchors are
**linearly interpolated**. Unless a rule below states otherwise, missing
applicable evidence is treated as **neutral rather than negative**,
usually receiving 50% of the available points.

Stage 6 evaluates transcript- and gene-specific hypotheses separately.
Evidence belonging to different genes is not combined. For each gene,
the strongest transcript-level pathogenicity hypothesis is retained, and
the final selected hypothesis determines the reported scored gene,
transcript, route, component scores, and final score.

------------------------------------------------------------------------

# 1. Variant Pathogenicity Score

## 1.1 Overview

The **Variant Pathogenicity Score** ranges from **0 to 100**.

For nuclear variants it consists of:

  Evidence group                          Maximum points
  ------------------------------------- ----------------
  CADD                                                20
  Conservation                                        10
  Route-specific evidence                             70
  **Total before ClinVar adjustment**            **100**

A variant may support more than one pathogenicity route. Each applicable
route is scored independently for each transcript, and the
**highest-scoring route for that transcript** is retained. Within each
gene, the transcript with the highest complete pathogenicity score is
used.

After the 0--100 pathogenicity score is calculated, a **ClinVar
adjustment** is applied. The final pathogenicity score is then bounded
to the range **0--100**.

### Missing evidence

Missing predictor values generally receive a neutral score equal to half
of the predictor's available contribution. Missing evidence is therefore
not interpreted as benign.

The accompanying **Pathogenicity Evidence Coverage** records how much of
the applicable pathogenicity evidence was actually observed. A neutral
score assigned because evidence is missing does not count as observed
evidence.

------------------------------------------------------------------------

## 1.2 CADD --- 0 to 20 points

CADD PHRED is converted first to a **0--10 evidence score**, then
multiplied by two.

    CADD PHRED   Evidence score   Pathogenicity points
  ------------ ---------------- ----------------------
            ≤5                0                      0
            10                2                      4
            15                4                      8
            20                6                     12
          22.5                8                     16
           ≥25               10                     20

Values between anchors are linearly interpolated.

**Missing CADD:** 5/10 → **10/20 points**.

------------------------------------------------------------------------

## 1.3 Conservation --- 0 to 10 points

Conservation is assessed using **phyloP100way** and **phastCons100way**.

When both are available, their 0--10 evidence scores are averaged. When
only one is available, the available score is used. When neither is
available, conservation receives the neutral score of **5/10**.

### phyloP100way scale

    phyloP100way   Score
  -------------- -------
            ≤0.0       0
             0.5       2
             1.0       4
             1.5       6
             2.0       8
            ≥2.5      10

Intermediate values are linearly interpolated.

### phastCons100way scale

    phastCons100way   Score
  ----------------- -------
              ≤0.00       0
               0.15       2
               0.30       4
               0.45       6
               0.60       8
              ≥0.75      10

Intermediate values are linearly interpolated.

**Both missing:** **5/10 points**.

------------------------------------------------------------------------

# 1.4 Missense Route --- 0 to 70 points

The missense route combines:

  Evidence                   Maximum points
  ------------------------ ----------------
  REVEL                                  30
  AlphaMissense                          30
  Protein/domain context                 10
  **Total**                          **70**

## 1.4.1 REVEL --- 0 to 30 points

REVEL is converted to a 0--10 score and multiplied by three.

    REVEL   Evidence score   Route points
  ------- ---------------- --------------
    ≤0.18                0              0
     0.30                2              6
     0.50                4             12
     0.64                6             18
     0.77                8             24
    ≥0.93               10             30

Intermediate values are linearly interpolated.

**Missing REVEL:** 5/10 → **15/30 points**.

## 1.4.2 AlphaMissense --- 0 to 30 points

AlphaMissense is converted to a 0--10 score and multiplied by three.

    AlphaMissense   Evidence score   Route points
  --------------- ---------------- --------------
            ≤0.15                0              0
             0.30                4             12
             0.45                6             18
             0.60                8             24
            ≥0.75               10             30

Intermediate values are linearly interpolated.

**Missing AlphaMissense:** 5/10 → **15/30 points**.

## 1.4.3 Protein/domain context --- 0 to 10 points

  ---------------------------------------------------------------------
  Protein context                                                Points
  ----------------------------- ---------------------------------------
  Annotated protein domain is                                        10
  present                       

  Amino-acid/protein-position                                         5
  information is present, but   
  no domain is identified       

  Protein context is present                                          0
  but provides neither domain   
  nor protein-position support  

  Missing                                                             5
  ---------------------------------------------------------------------

------------------------------------------------------------------------

# 1.5 Splice Route --- 0 to 70 points

The splice route combines:

  Evidence                                Maximum points
  ------------------------------------- ----------------
  SpliceAI                                            40
  Transcript-specific splice distance                 30
  **Total**                                       **70**

## 1.5.1 SpliceAI --- 0 to 40 points

SpliceAI is converted to a 0--10 score and multiplied by four.

    SpliceAI   Evidence score   Route points
  ---------- ---------------- --------------
        0.00                0              0
        0.10                2              8
        0.20                4             16
        0.50                6             24
        0.65                8             32
       ≥0.80               10             40

Intermediate values are linearly interpolated.

**Missing SpliceAI:** 5/10 → **20/40 points**.

## 1.5.2 Splice-junction distance --- 0 to 30 points

Absolute distance from the relevant splice boundary is used.

    Distance   Evidence score   Route points
  ---------- ---------------- --------------
       ≤5 bp               10             30
    6--20 bp                5             15
     \>20 bp                0              0
     Missing                5             15

------------------------------------------------------------------------

# 1.6 Loss-of-Function Route --- 0 to 70 points

The LoF route combines:

  Evidence                      Maximum points
  --------------------------- ----------------
  LoF / NMD context                         40
  Gene LOEUF                                20
  MANE transcript relevance                 10
  **Total**                             **70**

## 1.6.1 LoF / NMD context --- 0 to 40 points

  Evidence                                              Evidence score   Route points
  --------------------------------------------------- ---------------- --------------
  LOFTEE high-confidence or predicted/triggered NMD                 10             40
  LOFTEE low-confidence or NMD escape                                5             20
  Explicit evidence that the event is not LoF                        0              0
  Other available but non-decisive LoF/NMD evidence                  5             20
  Both LoF and NMD evidence missing                                  5             20

## 1.6.2 LOEUF --- 0 to 20 points

Lower LOEUF indicates stronger intolerance to loss-of-function
variation.

         LOEUF   Evidence score   Route points
  ------------ ---------------- --------------
        \<0.35               10             20
    0.35--0.90                5             10
        \>0.90                0              0
       Missing                5             10

## 1.6.3 MANE transcript relevance --- 0 to 10 points

  ---------------------------------------------------------------------
  Transcript status                                              Points
  ----------------------------- ---------------------------------------
  Scored transcript is a listed                                      10
  MANE transcript               

  A MANE transcript is                                                5
  available for the variant,    
  but the scored transcript is  
  not MANE                      

  No MANE transcript                                                  5
  information is available      
  ---------------------------------------------------------------------

------------------------------------------------------------------------

# 1.7 In-Frame Indel Route --- 0 to 70 points

The in-frame route combines:

  Evidence                                Maximum points
  ------------------------------------- ----------------
  Functional-site / in-frame evidence                 50
  Protein/domain context                              20
  **Total**                                       **70**

## 1.7.1 Functional-site evidence --- 0 to 50 points

  ---------------------------------------------------------------------
  Evidence                      Evidence score             Route points
  ------------------- ------------------------ ------------------------
  Active site,                              10                       50
  binding site,                                
  functional site, or                          
  explicit                                     
  deleterious support                          

  Explicit no-support                        0                        0
  / not-supporting /                           
  benign evidence                              

  Other available                            5                       25
  non-decisive                                 
  evidence                                     

  Missing                                    5                       25
  ---------------------------------------------------------------------

## 1.7.2 Protein/domain context --- 0 to 20 points

The same protein-context scale used for missense variants is applied,
multiplied by two.

  Protein context                                          Route points
  ------------------------------------------------------ --------------
  Annotated protein domain                                           20
  Protein/amino-acid position without domain                         10
  Available context without domain or position support                0
  Missing                                                            10

------------------------------------------------------------------------

# 1.8 Regulatory / Noncoding Route --- 0 to 70 points

The regulatory route combines:

  Evidence                       Maximum points
  ---------------------------- ----------------
  Regulatory annotation                      30
  Target-gene / TSS distance                 20
  ReMM                                       20
  **Total**                              **70**

## 1.8.1 Regulatory annotation --- 0 to 30 points

  -----------------------------------------------------------------------
  Regulatory evidence             Evidence score             Route points
  --------------------- ------------------------ ------------------------
  Promoter / PLS                              10                       30
  annotation                                     

  pELS, dELS, or                               5                       15
  enhancer annotation                            

  Explicitly not                               0                        0
  regulatory-relevant                            

  Regulatory-relevant                          5                       15
  but without a                                  
  stronger recognized                            
  annotation                                     

  Missing regulatory                           5                       15
  context                                        
  -----------------------------------------------------------------------

## 1.8.2 TSS distance --- 0 to 20 points

Absolute transcript-specific distance from the transcription start site
is used.

    TSS distance   Evidence score   Route points
  -------------- ---------------- --------------
           ≤1 kb               10             20
       \>1--2 kb                8             16
       \>2--5 kb                6             12
      \>5--10 kb                4              8
         \>10 kb                2              4
         Missing                5             10

## 1.8.3 ReMM --- 0 to 20 points

            ReMM   Evidence score   Route points
  -------------- ---------------- --------------
          \<0.50                0              0
    0.50--\<0.75                2              4
    0.75--\<0.90                4              8
    0.90--\<0.95                6             12
    0.95--\<0.98                8             16
           ≥0.98               10             20
         Missing                5             10

------------------------------------------------------------------------

# 1.9 ClinVar Pathogenicity Adjustment

ClinVar modifies the completed variant pathogenicity score additively.
Only a main classification consisting of Pathogenic, Likely pathogenic, or both
qualifies for P/LP rescue, scoring bonuses, priority and the P_LP report tab.
Equivalent spelling and separators are accepted. Conflicts, mixed non-P/LP
labels, missing values and descriptive text do not qualify; source evidence is
preserved. Gene-level ClinVar summaries are evaluated separately.

  ---------------------------------------------------------------------
  ClinVar evidence                                           Adjustment
  ----------------------------- ---------------------------------------
  Pathogenic only                                                   +20

  Pathogenic + Likely                                               +15
  Pathogenic                    

  Pathogenic mixed with                                               0
  conflicting, other, VUS,      
  benign, or likely benign      
  evidence                      

  Likely Pathogenic only                                            +15

  Likely Pathogenic mixed with                                        0
  conflicting, other, VUS,      
  benign, or likely benign      
  evidence                      

  VUS without conflicting or                                         +5
  other evidence                

  VUS with conflicting or other                                       0
  evidence                      

  Benign only                                                       −20

  Benign + Likely Benign                                            −15

  Benign mixed with                                                 −10
  conflicting, other, or VUS    
  evidence                      

  Likely Benign only                                                −15

  Likely Benign mixed with                                          −10
  conflicting, other, or VUS    
  evidence                      

  Absent or other unsupported                                         0
  classification                
  ---------------------------------------------------------------------

After this adjustment, the Variant Pathogenicity Score is restricted to
**0--100**.

------------------------------------------------------------------------

# 1.10 Pathogenicity Evidence Coverage

Pathogenicity Evidence Coverage ranges from **0 to 1** and represents
the fraction of the applicable pathogenicity evidence for which an
observed value was available.

For the common evidence layer:

-   observed CADD contributes 20% coverage;
-   conservation contributes 10% coverage when at least one conservation
    source is available.

The applicable route contributes its corresponding 70% according to
which route-specific inputs were observed.

Neutral points assigned to missing evidence affect the pathogenicity
score but **do not increase evidence coverage**.

------------------------------------------------------------------------

# 1.11 Mitochondrial Pathogenicity

A mitochondrial-specific pathogenicity model is **not currently
implemented in Stage 6**. Mitochondrial candidates continue with partial
scoring using available COMMON evidence (CADD, conservation and the main
ClinVar adjustment), plus existing inheritance, gene-disease and case
relevance rules. Nuclear route-specific predictors do not substitute for
an MT model. Missing evidence follows the existing defaults; weights are
unchanged. The scored route is COMMON, while source MT labels are retained.

MT candidates without transcript annotation can still receive partial
evaluation; the transcript stays empty. A single known source gene can
provide gene/case evidence. Missing gene context uses existing zero defaults.
Malformed identities and contradictory transcript references still fail QC.

MT rarity is unavailable for scoring. MT frequencies are preserved for
review and are not interpreted as nuclear carrier counts. The existing
missing AF/count defaults give `Rarity_score=50`, contributing five points
to the final score. This is a missing-evidence default, not established
rarity support. The workbook marks partial MT evaluation and explains the
default in score-cell comments. Existing worksheet thresholds still apply;
all retained candidates remain in ALL.

------------------------------------------------------------------------

# 2. Mode of Inheritance Score

## 2.1 Overview

The **Inheritance Score** ranges from **0 to 100** and measures
compatibility between the candidate's inheritance model and the
inheritance model or models expected for the case.

  ---------------------------------------------------------------------
  Relationship to expected                                        Score
  inheritance                   
  ----------------------------- ---------------------------------------
  Direct match                                                      100

  Matching AR                                                        80
  compound-heterozygous model   

  No expected inheritance                                            50
  information                   

  Candidate inheritance                                              50
  information unavailable       

  AD/AD de novo alternative                                          25
  from a compound candidate     
  when an AR model is expected  

  Incompatible                                                        0
  ---------------------------------------------------------------------

A direct compound-heterozygous match receives **80 rather than 100**.

Inheritance aliases are interpreted consistently, including X-linked
hemizygous as X-linked and mtDNA as mitochondrial.

------------------------------------------------------------------------

# 3. Gene / Disease Evidence Score

## 3.1 Overview

The **Gene / Disease Evidence Score** ranges from **0 to 100** and
combines four evidence groups:

  Evidence group                         Maximum points
  ------------------------------------ ----------------
  Gene--disease association strength                 50
  Mode-of-inheritance compatibility                  15
  ClinVar gene landscape                             20
  Literature support                                 15
  **Total**                                     **100**

If no gene can be assigned to the hypothesis, the Gene / Disease
Evidence Score is **0**.

------------------------------------------------------------------------

## 3.2 Gene--Disease Association Strength --- 0 to 50 points

The strongest recognized disease-association evidence is converted as
follows:

  Association classification        Points
  ------------------------------- --------
  Definitive                            50
  Strong                                50
  Moderate                              25
  Supportive                            25
  Limited                               25
  Animal Model Only                     25
  Disputed                               0
  Refuted                                0
  No Known Disease Relationship          0
  Missing or unrecognized               25

Missing or unrecognized evidence is neutral.

------------------------------------------------------------------------

## 3.3 Gene--Disease MOI Compatibility --- 0 to 15 points

The candidate inheritance model is compared with inheritance evidence
reported for the gene--disease relationship.

Candidate models are grouped into broad inheritance families:

-   AD, AD de novo, and Mosaic → autosomal-dominant family;
-   AR homozygous and AR compound heterozygous → autosomal-recessive
    family;
-   X-linked → X-linked family;
-   Mitochondrial → mitochondrial/maternal family.

  ---------------------------------------------------------------------
  Compatibility                                                  Points
  ----------------------------- ---------------------------------------
  Matching inheritance family                                        15

  Mosaic candidate with                                             7.5
  reported AD inheritance       

  Candidate family or                                               7.5
  gene--disease inheritance     
  evidence unavailable          

  Incompatible inheritance                                            0
  family                        
  ---------------------------------------------------------------------

------------------------------------------------------------------------

# 3.4 ClinVar Gene Landscape --- 0 to 20 points

The ClinVar gene landscape has three parts:

  Evidence                                          Maximum points
  ----------------------------------------------- ----------------
  Matching pathogenic molecular class                           10
  Number of P/LP variants in the gene                            5
  Fraction of classified variants that are P/LP                  5

## 3.4.1 Matching pathogenic molecular class --- 0 to 10 points

The selected pathogenicity route is compared with the gene's known P/LP
molecular-consequence landscape.

Examples of route matching include:

-   Missense → missense;
-   Splice → splice;
-   LoF → frameshift, nonsense, start-loss, or stop-loss;
-   In-frame → in-frame;
-   Regulatory → regulatory, intronic, UTR, noncoding, or near-gene.

  Evidence                                                    Points
  --------------------------------------------------------- --------
  At least one P/LP variant of a matching molecular class         10
  No matching class, but the gene has other P/LP variants          5
  Landscape missing and P/LP count is missing or \>0               5
  No known P/LP variants supporting the gene                       0
  Route has no defined landscape class                             5

## 3.4.2 P/LP variant count --- 0 to 5 points

    Gene P/LP count   Points
  ----------------- --------
                  0        0
                  1        1
               2--4        2
              5--10      3.5
               \>10        5
            Missing      2.5

## 3.4.3 Pathogenic fraction --- 0 to 5 points

    P/LP fraction   Points
  --------------- --------
                0        0
      \>0--\<0.05        1
     0.05--\<0.10      2.5
            ≥0.10        5
          Missing      2.5

------------------------------------------------------------------------

# 3.5 Literature Support --- 0 to 15 points

Literature evidence is scored from the previously collected publication
summary.

Only a completed, validated search can establish explicitly zero reviewed
publications. If searches exhaust their retries, available partial evidence
is retained and labeled incomplete; unavailable evidence follows the
existing missing-evidence rule. Validation and retry handling do not change
the literature scoring weights.

  Literature evidence                                           Points
  ----------------------------------------------------------- --------
  ≥2 human-case publications                                        15
  1 human-case publication                                          12
  Functional or animal evidence, without human-case support        7.5
  Explicitly zero reviewed publications / no publications            0
  Literature information missing                                   7.5
  Available but not classifiable into the categories above         7.5

Human-case evidence therefore receives the strongest literature
contribution.

------------------------------------------------------------------------

# 4. Case Relevance Score

## 4.1 Overview

The **Case Relevance Score** ranges from **0 to 100** and measures how
strongly the scored gene is connected to the specific clinical case.

Evidence is evaluated hierarchically. A predefined case-gene match or
direct disease-name match receives full support immediately. Otherwise,
broad disease-domain and HPO evidence are evaluated, and the strongest
applicable support is retained.

If no gene is assigned to the hypothesis, the Case Relevance Score is
**0**.

------------------------------------------------------------------------

## 4.2 Predefined Case Gene --- 100 points

A gene already identified as a case gene from the supplied case context
receives:

**100 points**

This includes genes already designated through the case-gene information
supplied to the workflow.

------------------------------------------------------------------------

## 4.3 Direct Disease Match --- 100 points

If the gene's disease-association evidence directly contains a disease
name supplied for the case:

**100 points**

------------------------------------------------------------------------

## 4.4 Broad Disease-Domain Match --- 70 points

When broad disease domains are supplied and the gene's
disease-association evidence contains a matching domain:

**70 points**

This evidence is weaker than a direct disease-name match.

------------------------------------------------------------------------

## 4.5 HPO Overlap --- 0 to 100 points

HPO support is calculated only when both patient HPO terms and
gene-associated HPO terms are available.

The score is:

**100 × (number of matched gene HPO terms / number of gene-associated
HPO terms)**

Matched HPO terms are used when available; otherwise, direct overlap
between patient HPO IDs and gene HPO IDs is used.

Examples:

    Matched gene HPO terms   Total gene HPO terms   HPO score
  ------------------------ ---------------------- -----------
                         0                     10           0
                         2                     10          20
                         5                     10          50
                         8                     10          80
                        10                     10         100

When both HPO and broad-domain evidence are applicable, the **higher
score** is retained.

If patient HPO terms or gene HPO terms are unavailable, HPO does not
contribute a negative score.

------------------------------------------------------------------------

## 4.6 No Supporting Case Connection

If case-specific information exists but none of the available evidence
supports a connection between the gene and the case:

**20 points**

If **no case-specific information is available at all**:

**50 points**

The latter represents an unknown/neutral state rather than evidence
against the gene.

------------------------------------------------------------------------

# 5. Rarity Score

## 5.1 Overview

The **Rarity Score** ranges from **0 to 100** and ranks variants that
have already passed the earlier rarity filters.

It combines:

  Evidence                                             Weight within Rarity Score
  -------------------------------------------------- ----------------------------
  Population allele-frequency support                                         80%
  Inheritance-appropriate population count support                            20%

Thus:

**Rarity Score = 10 × (0.8 × AF evidence score + 0.2 × count evidence
score)**

Both underlying evidence scores use a 0--10 scale.

The effective population AF is the **higher** of the available overall
gnomAD AF and ancestry-group maximum AF.

------------------------------------------------------------------------

## 5.2 AF Thresholds by Inheritance Model

  Candidate model                             AF reference threshold
  ----------------------------------------- ------------------------
  AD                                                          1×10⁻⁴
  AD in a case gene                                           1×10⁻³
  AR compound heterozygous                                    1×10⁻³
  AR compound heterozygous in a case gene                     1×10⁻²
  AR homozygous                                               1×10⁻³
  X-linked                                                    5×10⁻³
  Mosaic                                                      1×10⁻³
  Mitochondrial                                               5×10⁻³

------------------------------------------------------------------------

## 5.3 AF Rarity Scale --- 0 to 10

AF support is determined by the ratio between the applicable threshold
and the observed effective AF. The scale is logarithmic, with linear
interpolation between the log-scaled anchors.

  Observed AF relative to threshold     AF evidence score
  ----------------------------------- -------------------
  10× above threshold                                   0
  Exactly at threshold                                  1
  5× below threshold                                    5
  50× below threshold                                   8
  ≥500× below threshold                                10
  Confirmed AF = 0                                     10
  Missing AF                                            5

For example, under an AD threshold of 1×10⁻⁴:

    Effective AF Relationship            AF score
  -------------- --------------------- ----------
          1×10⁻³ 10× above threshold            0
          1×10⁻⁴ threshold                      1
          2×10⁻⁵ 5× below                       5
          2×10⁻⁶ 50× below                      8
         ≤2×10⁻⁷ ≥500× below                   10

If AF is missing but all relevant nuclear population genotype counts are
explicitly present and equal to zero, the nuclear variant is treated as
absent from the population and receives the maximum AF support.

------------------------------------------------------------------------

# 5.4 Population Count Support

The population count used depends on the inheritance model.

## 5.4.1 Heterozygous-count scale

Used for AD, AD de novo, Mosaic, and as part of compound-heterozygous
scoring.

    gnomAD heterozygous count   Score
  --------------------------- -------
                            0      10
                         1--5       8
                        6--15       6
                       16--30       4
                       31--49       2
                          ≥50       0
                      Missing       5

## 5.4.2 Homozygous / hemizygous-count scale

Used for AR homozygous and X-linked hemizygous support, and as part of
compound-heterozygous scoring.

      Count   Score
  --------- -------
          0      10
          1       6
          2       2
         ≥3       0
    Missing       5

For X-linked candidates, the hemizygous population count is evaluated
with this scale.

## 5.4.3 Compound-heterozygous allele count

For each compound-heterozygous allele, count support is the mean of:

-   the heterozygous-count score; and
-   the homozygous-count score.

Each allele therefore receives its own rarity score.

------------------------------------------------------------------------

# 6. Compound-Heterozygous Scoring

## 6.1 Allele-Level Scoring

Each member of a valid compound-heterozygous pair is scored as a
complete candidate using:

-   its own Variant Pathogenicity Score;
-   compound-heterozygous Inheritance Score;
-   gene/disease evidence for the shared gene;
-   case relevance for the shared gene;
-   its own compound-allele Rarity Score.

A matching compound-heterozygous inheritance model receives an
**Inheritance Score of 80**.

------------------------------------------------------------------------

## 6.2 Compound Pair Score

For a selected pair:

**Compound Pair Score = mean of the two alleles' Final Scores**

The pair score is therefore also on a **0--100 scale**.

Only the selected pair receives a Compound Pair Score and an explicit
selected partner. Additional retained alleles scored against the same AR
compound-candidate gene retain AR inheritance, gene-disease compatibility
and allele-specific rarity even when they are outside the selected main
pair. AD fallback applies only when an AD hypothesis for a different gene
wins for that variant. Stage 6 compares fully scored gene hypotheses before
selection; weak or missing evidence uses the existing component rules.
Source compound relationships are retained, while selected partner and pair
score are blank for AD fallbacks. Stage 7 displays the selected inheritance;
the AR_compound tab still requires a reciprocal selected pair scoring at
least 30.

------------------------------------------------------------------------

# 7. Final Score and Hypothesis Selection

## 7.1 Final Score

For every complete candidate hypothesis:

**Final Score = 0.30 × Variant Pathogenicity Score\
+ 0.20 × Inheritance Score\
+ 0.25 × Gene / Disease Evidence Score\
+ 0.15 × Case Relevance Score\
+ 0.10 × Rarity Score**

The result is reported on a **0--100 scale** and rounded to two decimal
places.

------------------------------------------------------------------------

## 7.2 Transcript and Route Selection

For each transcript, every applicable pathogenicity route is evaluated
independently.

The route with the highest Variant Pathogenicity Score becomes the
**Scored Route** for that transcript.

For each gene, the transcript with the highest Variant Pathogenicity
Score becomes the gene's selected transcript hypothesis.

This preserves transcript-specific scoring while preventing evidence
from different transcripts or genes from being combined into a single
artificial hypothesis.

------------------------------------------------------------------------

## 7.3 Gene-Hypothesis Selection

When a variant is associated with more than one gene, each gene receives
an independent complete scoring hypothesis.

Case-gene hypotheses and hypotheses with a pure Pathogenic, Likely
pathogenic or combined P/LP main classification receive selection priority. Within
the same priority class, hypotheses are compared by:

1.  Final Score;
2.  Variant Pathogenicity Score;
3.  Gene / Disease Evidence Score.

The selected hypothesis determines the reported scored gene, transcript,
route, component scores, and Final Score.

------------------------------------------------------------------------

## 7.4 Compound Pair Selection

Valid compound pairs are compared with the alternative single-variant
hypotheses for their two alleles.

A compound hypothesis is retained when it has the applicable
case-gene/P/LP priority or, within the same priority class, when its
pair score is at least as high as the competing alternative-gene score for
both alleles. An AR allele hypothesis in the same gene supports that pair;
its individual score does not compete against the pair's mean score.

Only one compound pair is selected per gene.

### Fallback to a Stronger Alternative Gene Hypothesis

A variant that participates in a compound-heterozygous pair may also
support an independent hypothesis in another gene. The compound
assignment does not automatically override that alternative hypothesis.

If either allele has a stronger competing gene hypothesis, according to
the same hypothesis-selection priority and score comparison rules, the
compound pair is not selected over that stronger alternative. The
affected variant therefore falls back to its stronger independently
scored gene hypothesis rather than being forced into the compound pair.

For example, if variants A and B form a compound pair in gene X with a
Compound Pair Score of 72, but variant A has an alternative hypothesis
in gene Y with a score of 81 within the same priority class, the
compound pair in gene X does not replace the stronger gene Y hypothesis
for variant A.

This comparison is performed for both members of the pair: a
non-priority compound pair must be at least as strong as the competing
singleton hypothesis of each allele in order to be selected.

------------------------------------------------------------------------

# 8. Missing Evidence

The scoring system distinguishes **missing evidence** from **negative
evidence**.

As a general rule, when an applicable scoring input is missing, it
receives a **neutral value equal to 50% of the points available for that
input**. Missing evidence therefore does not behave as evidence against
the candidate.

The main exceptions are:

-   **Conservation:** if both phyloP and phastCons are missing, the
    conservation component receives the neutral value of 5/10. If only
    one of the two is available, the available conservation score is
    used rather than averaging it with a neutral value.
-   **Regulatory annotation/context:** missing regulatory context is
    handled by the route-specific neutral rule defined in the Regulatory
    Pathogenicity section.
-   **Case relevance:** 50/100 is used only when no case-specific
    information is available. If case information exists but does not
    support a connection between the candidate gene and the case, the
    score can be lower.
-   **Population AF:** if AF itself is missing but the relevant nuclear
    population genotype counts are explicitly available and all are
    zero, the variant can be treated as absent from the population
    dataset rather than receiving the neutral missing-AF score.

Neutral points assigned because evidence is missing should not be
interpreted as positive biological evidence. Explicit non-supporting or
contradictory evidence can reduce the corresponding component below the
neutral value, including to zero where defined by the scoring scale.

------------------------------------------------------------------------

# 9. Stage 6 Scoring Outputs

Stage 6 reports the following scoring results for the selected
hypothesis:

-   **Variant Pathogenicity Score** --- 0--100 pathogenicity score for
    the selected transcript and route.
-   **Pathogenicity Evidence Coverage** --- fraction from 0--1
    indicating how much applicable pathogenicity evidence was observed.
-   **Inheritance Score** --- 0--100 compatibility with expected case
    inheritance.
-   **Gene / Disease Evidence Score** --- 0--100 gene/disease support.
-   **Case Relevance Score** --- 0--100 case-specific relevance.
-   **Rarity Score** --- 0--100 population rarity support.
-   **Scored Gene** --- gene hypothesis used for the retained score.
-   **Scored Transcript** --- transcript used for the retained
    pathogenicity score.
-   **Scored Route** --- pathogenicity route selected for that
    transcript.
-   **Scored Inheritance** --- inheritance model actually used for scoring.
-   **AD Fallback** --- identifies a compound candidate recalculated as AD.
-   **Scored Compound Partner** --- selected partner when the variant
    belongs to the chosen compound pair.
-   **Compound Pair Score** --- mean Final Score of the two selected
    compound alleles.
-   **Final Score** --- weighted 0--100 prioritization score.

The scores are intended to determine **review priority**. They summarize
the evidence available to the workflow and do not replace manual
biological or clinical interpretation.

For Stage 7 **Top_candidates**, every non-case-gene variant must have
Variant Pathogenicity Score ≥ 30 in addition to the existing final-score
and expected-inheritance requirements. This minimum applies separately
to each compound allele; a strong partner cannot waive it. Case-gene
variants retain their existing eligibility without this additional minimum.
