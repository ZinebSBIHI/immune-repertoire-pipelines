# Design notes

Why this repository makes the choices it does. Each section states a
decision, the failure it avoids, and where in the code it lives. These are
the places where a repertoire analysis can be technically correct and
scientifically wrong, so they are documented rather than buried in a default.

Nothing here is a rule about what a study must conclude. They are rules about
what an analysis is entitled to claim.

---

## 1. Sequencing depth confounds almost everything

**Decision.** Richness, clone counts and clonal concentration are compared
only at a common depth, obtained by resampling.
*`rarefy_to_common_depth()` in `R/rep-rarefaction.R`.*

Observed richness is a monotone function of depth: sequence deeper and you
find more clonotypes, with no biological change whatever. The same holds for
the number of expanded clones, and inversely for the fraction of the
repertoire held by the largest clone. When sequencing depth differs between
the groups being compared — and it usually does, because depth depends on
input material, which depends on the condition — a comparison of these
quantities measures the depth difference.

Two choices inside this:

- **The common depth is taken within subject, not across the cohort.** For a
  within-subject comparison of arms or timepoints, degrading every subject
  to the depth of the shallowest sample in another subject discards
  information for nothing. `within = "subject_id"` is the default.
- **Resampling is without replacement.** The question is what would have
  been seen with fewer molecules from the same library, which is sampling
  from a finite pool, not from an estimated distribution.

Metrics differ in how exposed they are. Morisita-Horn overlap is
abundance-weighted and comparatively robust; Jaccard counts clonotypes and
is not. Where depths differ and only one can be reported, report
Morisita-Horn.

---

## 2. The abundance floor is a decision, not a default

**Decision.** `repertoire$min_abundance` is swept before it is fixed, and
the cost of each value is reported.
*`sweep_abundance_floor()` in `R/ana-differential.R`; section 02 of the
bulk pipeline.*

The floor below which a clonotype is not counted is the most consequential
free parameter in bulk repertoire work, and it acts on de novo calls in two
opposing directions at once:

- raising it removes single-molecule clonotypes, which reduces the calls;
- raising it also pushes low-abundance *baseline* clonotypes below the
  detection threshold, so they count as "not detected at baseline", which
  **creates** de novo calls.

Which effect dominates is a property of the dataset, not of the method. On
the simulated example shipped with this repository the second dominates, and
the number of de novo calls *rises* by about 30% when the floor moves from
1 to 3, while a quarter of the clonotypes are discarded. An analysis that picks a floor by habit and reports the result
has no way of knowing which way its choice pushed the answer.

Hence: run the sweep, fix the value, state it in Methods, and include the
sweep as a supplementary figure. The default here is 1, the most inclusive
and least specific setting, so that the filtering is visible as something
added rather than something assumed.

A related point of hygiene: UMI counts and read counts are not
interchangeable, and `duplicate_count` in an AIRR file means different
things depending on the assembler. The abundance column is therefore named
explicitly and never inferred. *`collapse_clonotypes()` raises an error
rather than guessing.*

---

## 3. Category and significance are different questions

**Decision.** Differential abundance returns two independent axes and never
merges them into a single causal label.
*`test_clonal_abundance()` in `R/ana-differential.R`.*

- **Category** is descriptive, relative to the subject's own baseline:
  `de_novo`, `pre_existing_expanded`, `stable`, `contracting`, `lost`.
- **Significance** is statistical: a Fisher exact test on the clonotype's
  molecule count against the rest of the repertoire, at both timepoints,
  with FDR control.

The two do not coincide. A clonotype can be de novo and non-significant —
one molecule, absent before, where the test has no power and where absence
at baseline is exactly what a detection floor produces by chance. It can be
pre-existing and strongly significant. The pipeline reports a combined
`status` that names both axes side by side, and nothing that asserts
induction.

The categorical scheme follows Sahin et al. (*Nature* 2026,
[doi:10.1038/s41586-025-10004-2](https://doi.org/10.1038/s41586-025-10004-2),
PMID 41708868). What that study did *not* do is as important as what it did:
candidate selection there combined bulk tracking with paired single-cell
TRA/TRB sequencing and functional validation, and used only the two
categories. Everything beyond that in this implementation is this
repository's own, and the word "induced" is deliberately absent from the
output. Calling a clonotype vaccine-induced on the strength of a bulk
fold-change claims the functional validation that was not performed.

---

## 4. Convergence needs a null

**Decision.** The number of public clonotype groups is always reported
against a label-permutation null.
*`permute_public_groups()` in `R/ana-convergence.R`.*

The count of public groups rises with repertoire depth, with the number of
subjects, and with a looser identity threshold. On its own it is not
interpretable. Permuting subject labels while holding the sequences fixed
preserves sequence composition and V gene structure and randomises only the
assignment to subjects, which gives the number of groups expected with no
biological sharing.

This matters more than it may appear. Soto et al. (*Nature* 2019,
[doi:10.1038/s41586-019-0934-8](https://doi.org/10.1038/s41586-019-0934-8),
PMID 30760926) found 1-6% of heavy-chain clonotypes shared between any two
healthy individuals, and shared clonotypes present in cord blood as well as
adult repertoires. Sharing is the baseline condition, not the finding. On
the simulated bulk dataset, where six CDR3 families were planted in every
subject and everything else is random, the pipeline finds 22 public groups
at an identity threshold of 0.75, against 21.3 expected by chance (p = 0.47
against the null). The six planted families are the groups spanning all
three subjects; the sixteen others, each shared by two subjects, are what
chance alone produces at that threshold and depth. Without the permutation
all 22 would read as convergence.

Three further choices:

- **Distance is normalised by length.** An absolute edit distance of 1 is a
  different biological statement on a 6-residue CDR3 than on a 22-residue
  one, and on short CDR3s an absolute threshold often yields no edges at
  all, so the analysis silently reports no convergence.
- **Single linkage, with the group's mean identity reported.** A lineage of
  related CDR3s forms a chain, and requiring every pair in a group to meet
  the threshold would fragment it. The cost is that long chains can link
  distant sequences, so `mean_identity` and `min_identity` are returned per
  group and should be inspected.
- **Germline-proximal sequences are excluded by default (BCR).** Unmutated,
  unswitched sequences sit close to germline and therefore close to each
  other. A convergence analysis that includes naive sequences measures
  shared germline, not shared selection.
  *`filter_informative_sequences()`.*

---

## 5. Small cohorts have a p-value floor

**Decision.** Every paired test reports the smallest p-value its sample size
can attain, and direction concordance across subjects is reported as a
first-class result.
*`attainable_p_paired()`, `paired_comparison()` in `R/ana-statistics.R`.*

A paired signed-rank test on *n* pairs cannot return a two-sided p below
2·2⁻ⁿ, whatever the effect size:

| n | smallest attainable two-sided p |
|---|---|
| 3 | 0.25 |
| 4 | 0.125 |
| 5 | 0.0625 |
| 6 | 0.03125 |

Below n = 6 no result can reach p < 0.05. Reporting "not significant" in
that regime describes the test rather than the data; reporting an
uncorrected p that happens to fall under 0.05 after several metrics were
tried describes neither.

What carries evidence at this size is consistency of direction. A
difference in the same direction in every subject is a statement about
reproducibility that does not depend on a test the design cannot support.
`format_paired_result()` writes the sentence the result actually supports,
including the bound, rather than a bare significance verdict.

---

## 6. Gene usage at clonotype level and at molecule level

**Decision.** V and J usage is computed both ways and both are reported.
*`gene_usage(level = ...)` in `R/rep-clonotypes.R`.*

At molecule level — or at cell level in single-cell data — one large clone
dominates its sample's usage profile. An apparent V gene bias can be a
single expanded lineage, and a transcriptome-level measurement cannot
distinguish the two because a per-cell measurement also reflects a large
clone. Computing usage once per clonotype separates repertoire-wide bias
from clonal expansion. Where the two levels disagree, the clone is the
explanation, and the finding is about that clone rather than about the
repertoire.

---

## 7. One heavy chain per cell, and what to do with the rest

**Decision.** Barcodes with two well-supported heavy chains are dropped as
doublets rather than resolved; elsewhere the dominant chain represents the
cell, with class-switched isotypes preferred.
*`resolve_one_heavy_per_cell()` in `R/bcr-cells.R`.*

A 10x V(D)J library reports more than one heavy chain per barcode for
several different reasons, which deserve different treatment: ambient mRNA
and index hopping give a dominant contig and a minor one; a true doublet
gives two well-supported contigs; and IGHM and IGHD are co-transcribed in
naive and marginal-zone B cells. Keeping the highest-UMI contig without
distinguishing these converts doublets into cells and makes the isotype call
a function of transcript abundance.

**The isotype consequence is easy to miss.** Because one isotype is retained
per cell, an IgM-versus-IgD split across naive or marginal-zone cells
describes which transcript happened to dominate, not two populations. The
quantity that means something is the combined unswitched fraction, and
`summarise_isotypes()` reports it that way by default.

---

## 8. A dominant clone is profiled, not silently excluded

**Decision.** Clones above 1% of a subject's cells are profiled with the
criteria needed to judge them; nothing is excluded by the code.
*`profile_dominant_clones()` in `R/bcr-cells.R`.*

A single lineage can account for most of a compartment, and then every
repertoire-level statistic on that compartment follows it: gene usage, CDR3
length, mutation load. Two criteria separate real expansion from a technical
artefact:

- UMI support at the level of the subject's other cells — a clone built from
  low-UMI contigs is suspect;
- a paired light chain across its cells — a heavy-chain-only expansion is
  usually an assembly artefact.

**Intraclonal diversity is deliberately not a criterion.** A genuine
expansion over a few days produces identical daughter cells, so low
diversity within a clone is what real expansion looks like. Using diversity
as a filter removes the clearest cases.

Exclusion is a judgement to be made and stated, so the function returns a
`verdict` recording which criteria each clone meets and leaves the decision
to the analyst. Any statistic that an exclusion changes should be reported
both ways, as a sensitivity analysis.

A related threshold: expansion cut-offs differ between resting and
stimulated material. A multi-day stimulation culture produces identical
daughter cells, so the cut-off that marks genuine expansion in resting
material is met trivially after stimulation. *`expansion_threshold()`.*

---

## 9. QC thresholds are per sample

**Decision.** Single-cell QC thresholds come from a per-sample block in the
configuration, and the report names every library that fell back on the
default.
*`section_gex_qc()` in `R/gex-core.R`; `gex$qc` in the configuration.*

A resting library and a stimulated one do not share an upper bound on
`nFeature_RNA`: proliferating cells are larger and carry more transcripts. A single global threshold therefore preferentially discards the
activated cells — which are usually the cells under study — while looking
like careful quality control.

---

## 10. Batch correction is off by default

**Decision.** `gex$integration$method` defaults to `none`, and the report
states whether correction was applied and on what variable.
*`section_gex_reduction()` in `R/gex-core.R`.*

Integration makes samples resemble each other. Where batch and condition
coincide — as they do whenever each condition is prepared as its own
library — that is indistinguishable from removing the effect under study.
Correction should be a decision taken because a batch effect was
demonstrated, and recorded as such.

With correction off, structure that follows subject or library stays visible
in the embedding instead of being removed from it, which is what makes the
next check possible.

---

## 11. A cluster dominated by one subject cannot support a group comparison

**Decision.** Every cluster's largest single-subject fraction is computed,
and clusters above the configured threshold are named in the report before
any proportion is compared.
*`section_cluster_composition()` in `R/gex-core.R`;
`gex$max_single_subject_fraction`.*

If a cluster is 95% one subject's cells, a difference in its size between
treatment arms is a difference between subjects. The comparison is not
rescued by a larger cell count: the cells are not independent observations
of the condition, they are one subject measured many times. Either exclude
such clusters from proportion comparisons and say so, or model the subject
effect explicitly.

---

## 12. Reproducibility is mechanical, not aspirational

Each pipeline run writes, into its own timestamped directory: the resolved
configuration, `sessionInfo()`, figure legends generated alongside the
figures, a table manifest, and a plain-text numeric report. Nothing
overwrites a previous run.

The reason is narrow and practical. A number quoted in a manuscript should
come from a file the analysis wrote, with the parameters and package
versions that produced it beside it — not from a console session that has
since been lost, and not from a figure directory whose contents were
replaced three revisions ago. Run directories that produced no output are
swept automatically, so what remains is what ran.

---

## 13. A database match is a hypothesis, and the matching rule is stated

**Decision.** Specificity annotation uses one stated rule per data type, and
an unannotated receptor is never read as evidence of anything.
*`R/ana-specificity-db.R`, `R/ana-specificity.R`, `R/sc-specificity.R`.*

Public receptor databases are dominated by a few intensively studied
pathogens — CMV, EBV, influenza and, since 2020, SARS-CoV-2 (Goncharov et
al. *Nat Methods* 2022, PMID 35970936). Most receptors specific for a less
studied antigen are simply absent from them, so a missing annotation says
nothing about a receptor's specificity, and a match is a sequence identity
with a published receptor, a hypothesis to test.

- **Bulk TCR: exact identity of the CDR3 amino-acid sequence**, against
  VDJdb and McPAS-TCR. The references are collapsed to one row per CDR3
  before joining, so a CDR3 recorded against several epitopes never
  duplicates a clonotype. In the species composition, each (CDR3, species)
  pair carries the clonotype's abundance, as in the original figure, so the
  fractions are of annotations rather than of molecules — stated in the
  figure legend.
- **Single-cell TCR: exact identity, then propagation within the dataset.**
  This is the rule of `Trex::annotateDB()`: a cell whose CDR3 is within *k*
  edits of a CDR3 that matched exactly *in the same data* takes that
  anchor's annotation, and the first anchor within reach in cell order wins,
  even when another is nearer. A reference CDR3 that no cell carries is
  never used. The rule is reimplemented rather than called, so that it runs
  on any reference and keeps every barcode intact; the test suite checks
  that it returns the same annotation as `annotateDB()` at *k* = 0, 1 and 2.
- **BCR: same IGHV gene and CDRH3 identity ≥ 0.80**, one minus the
  Levenshtein distance over the longer sequence, against a reference of
  antibodies with a known target such as CoV-AbDab (Raybould et al.
  *Bioinformatics* 2021, PMID 32805021). The conserved junction anchors are
  removed from the query, because the reference records CDRH3 without them,
  and V genes are reduced to `IGHV<family>-<number>` on both sides.

Among clonotypes that changed significantly, a match to a pathogen
unrelated to the study (a bystander) points to a change that is not driven
by the exposure studied. Categories are checked in a fixed order and the
first match wins; the antigen under study, when named, is checked first and
leaves the bystander list. The label is "significantly expanded", never
"induced": see §3.

---

## 14. Two views of change, kept apart

**Decision.** The detection-based view of new and expanded clonotypes is
reported beside the Fisher-based calls, never instead of them.
*`section_pharmacodynamics()` in `R/ana-pharmacodynamics.R`.*

The pharmacodynamic presentation used in clinical immunotherapy reports
(Weber et al. *Lancet* 2024, PMID 38246194) follows the focusing of the
repertoire and the appearance of clonotypes absent at baseline, with
detection thresholds (here at least 2 molecules and a frequency of 10⁻⁵)
rather than a test. It answers a different question from §3 — how much of
the repertoire is new, not which clonotypes changed beyond counting noise —
and the two are reported side by side.

Focusing is measured as the change from each subject's own baseline in
Shannon clonality and in the Gini index. **Shannon clonality normalises by
the observed richness, which grows with depth**, so where a subject's depth
varies between timepoints part of its trajectory is technical. A control
table gives the Spearman correlation of each metric with depth and the
agreement between the two metrics: a conclusion that holds for Gini as well
as for Shannon clonality does not depend on that choice.

---

## 15. The similarity threshold of a CDR3 network is a reported decision

**Decision.** CDR3 networks link sequences within a Levenshtein distance of
20% of the longer sequence, and the number of edges at neighbouring
thresholds is reported with them.
*`cdr3_similarity_network()`, `cdr3_network_diagnostic()` in `R/rep-cdr3.R`.*

An absolute distance of 1 is a far stronger statement about two 18-residue
CDR3s than about two 10-residue ones, and on real repertoires it often
yields no edge at all, while the edge count changes by orders of magnitude
between neighbouring fractions. The diagnostic table gives, for every
subject, the edges at 10%, 15% and 20% and at an absolute distance of 1, so
that the value used can be judged from the evidence rather than taken on
trust.

---

## 16. Following clonotypes across samples

**Decision.** Proportions are computed over each sample's whole repertoire
before any clonotype is selected, and the union of each sample's top
clonotypes is followed everywhere.
*`clone_compare_table()`, `section_clone_tracking()` in `R/rep-tracking.R`.*

Selecting first and normalising afterwards makes a clonotype's height
depend on which other clonotypes were selected and which samples are
compared. Over time, two sets are followed per subject: the 15 clonotypes
dominant at baseline, and the 15 significantly expanded clonotypes with the
largest fold change, separately for those new and those present at
baseline, with ties at the boundary kept. A clonotype not detected at a
timepoint is drawn at zero, never interpolated.

---

## 17. Selection and genealogy of B-cell clones

**Decision.** Selection is estimated on clonal consensus sequences, and
lineage trees are built only for the clones that link the sets of cells a
study compares.
*`R/bcr-selection.R`, `R/bcr-lineage.R`.*

- **Replacement-to-silent ratios** are reported for CDR and FWR separately:
  antigen selection favours replacement mutations in the CDRs and purges
  them from the framework (Hershberg et al. *Int Immunol* 2008,
  PMID 18397909).
- **BASELINe** (Yaari et al. *Nucleic Acids Res* 2012, PMID 22641856)
  collapses each clone to one sequence first, so that an expanded clone
  counts once and not once per cell, and uses the focused test statistic.
  Groups of fewer than five sequences or two clones are not tested, and are
  listed rather than dropped silently. The targeting model (HH_S5F) is
  named explicitly instead of left to a default.
- **Lineage trees** (maximum likelihood, rooted on the clone's germline;
  Hoehn et al. *PLoS Comput Biol* 2022) are built for clones present in
  every set of cells the configuration names — for example a sorted
  compartment and the cells that divided in culture — because those are the
  clones whose genealogy says how the sets are related.

---

## 18. Joint RNA and receptor neighbourhoods

**Decision.** The receptor embedding is combined with the transcriptome by
weighted nearest neighbours, and the weight each cell gives to its RNA is
reported.
*`section_tcr_embedding()` in `R/sc-embedding.R`.*

Trex turns each CDR3 into a vector (here with its geometric method, which
needs no deep-learning backend; Borcherding et al. *Nat Immunol* 2024,
PMID 39164479). Cells whose embedding is incomplete are removed before any
neighbour search. Weighted nearest neighbours (Hao et al. *Cell* 2021,
PMID 34062119) learn, cell by cell, how much of the neighbourhood the RNA
should define; the median RNA weight is reported, because a joint embedding
dominated by one modality is that modality's embedding. A CoNGA-style
reduction to one representative cell per clonotype (Schattgen et al.
*Nat Biotechnol* 2022, PMID 34426704) keeps expanded clonotypes from
dominating the receptor neighbourhoods.

---

## What this repository does not do

- **It does not decide what to exclude.** Dominant clones, donor-dominated
  clusters and low-coverage signatures are flagged; the judgement stays with
  the analyst, and is expected to be stated.
- **It does not infer specificity.** A CDR3 matching VDJdb or McPAS-TCR is a
  hypothesis. Database coverage is skewed towards a few well-studied
  epitopes, so the annotated fraction of a repertoire reflects what has been
  studied as much as what the repertoire recognises. An unannotated
  clonotype is uncharacterised, not non-specific.
- **It does not replace the upstream tools.** Cell Ranger, IgBLAST, MakeDb
  and MiXCR run outside R and are not reproduced; the pipelines start from
  their AIRR output.
- **It does not make a small cohort larger.** Where the design cannot reach
  a nominal significance level, the pipelines say so rather than reporting a
  p-value as though it could.
