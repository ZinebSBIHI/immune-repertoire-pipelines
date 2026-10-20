# immune-repertoire-pipelines

Reproducible R workflows for adaptive immune receptor repertoire (AIRR)
analysis and the single-cell transcriptome analyses that accompany it.

Four pipelines cover single-cell TCR repertoires, single-cell BCR
repertoires with clonal inference and somatic hypermutation, bulk TCR/BCR
repertoires with differential clonal abundance testing, and single-cell gene
expression. They are written to be applied to any cohort: sample labels,
thresholds and contrasts come from a single configuration file, and no path
or study-specific value is hard-coded in the analysis code.

What distinguishes this from a collection of analysis scripts is that the
decisions which make repertoire work go quietly wrong are made explicit and
diagnosable rather than left as defaults. Sequencing depth is equalised
before richness is compared. The abundance floor is swept and its cost
reported, because it pushes de novo calls in both directions at once.
Descriptive category and statistical significance are returned as separate
axes and never merged into a causal label. Convergence between subjects is
reported against a permutation null, because shared clonotypes are the
baseline condition rather than the finding. Paired tests report the smallest
p-value the sample size can attain. Each of these is argued, with the code
that implements it, in **[docs/design_notes.md](docs/design_notes.md)**.

---

## Quick start

The bulk pipeline runs end to end on simulated data, with no restricted
input, in a minute or two:

```bash
git clone https://github.com/ZinebSBIHI/immune-repertoire-pipelines.git
cd immune-repertoire-pipelines

# R packages for the bulk pipeline
Rscript -e 'install.packages(c("yaml","dplyr","tidyr","ggplot2","scales","stringdist"))'

# generate the example dataset
Rscript examples/simulate_airr.R

# configure and run
cp config/config.example.yml config/config.yml
Rscript pipelines/01_bulk_repertoire.R
```

The run writes a timestamped directory under `results/` holding figures,
tables, a plain-text numeric report, the resolved configuration and
`sessionInfo()`.

The simulated data carry no biological meaning. They exist so that the
pipeline can be read while it is running, and so that the methods can be
checked against a known answer: the simulator plants 25 expanded and 20
de novo clonotypes per subject, and the pipeline recovers them.

---

## Layout

The pipelines are named on a single axis — the kind of sequencing, which
fixes the input format and the unit of observation. The receptor (TCR or
BCR) is a configuration choice, not a separate pipeline, because the
repertoire metrics are identical; BCR simply adds the immunoglobulin-only
steps (clonal inference, germlines, somatic hypermutation, isotypes), which
switch on automatically for an Ig locus.

```
pipelines/                      one script per data type, driven by the config
  01_bulk_repertoire.R          bulk TCR or BCR   (unit: molecule / UMI)
  02_single_cell_repertoire.R   single-cell TCR or BCR   (unit: cell)
  03_single_cell_gex.R          single-cell transcriptome (Seurat v5)

R/
  setup.R                  locate the repository root, source the modules
  core/
    config.R               configuration loading and validation
    io.R                   run directories, figures, tables, reports
    plotting.R             theme and colour-blind-safe palettes
  repertoire/
    clonotypes.R           clonotype definition, collapsing, gene usage
    metrics.R              diversity, clonality, overlap
    rarefaction.R          resampling to a common depth
    sections.R             the shared analysis sections both repertoire
                           pipelines call, so they cannot drift apart
  analysis/
    differential_abundance.R  category and significance axes, floor sweep
    convergence.R             public clonotype groups and the permutation null
    statistics.R              paired tests at small n, attainable p
  bcr/
    cells.R                one heavy chain per cell, isotypes, dominant clones

config/config.example.yml  the only place study-specific values belong
docs/                      design notes and the verified bibliography
examples/simulate_airr.R   generates the example dataset
data/README.md             where to put inputs and reference databases
```

### Which pipeline for which data

| You have | Pipeline | Set in config |
|---|---|---|
| Bulk TCRβ (MiXCR / AIRR TSV, one file per sample) | `01_bulk_repertoire.R` | `locus: TRB` |
| Bulk BCR / IGH (AIRR TSV per sample) | `01_bulk_repertoire.R` | `locus: IGH` |
| Single-cell TCR (Cell Ranger VDJ contigs, or AIRR per library) | `02_single_cell_repertoire.R` | `locus: TRB` |
| Single-cell BCR (AIRR per library, MakeDb output) | `02_single_cell_repertoire.R` | `locus: IGH` |
| Single-cell gene expression (Cell Ranger `.h5`) | `03_single_cell_gex.R` | — |

The two repertoire pipelines run the same analysis sections — composition,
clonotype-definition sensitivity, structure, diversity at common depth, gene
usage, overlap, convergence against a permutation null, paired comparisons —
from one shared implementation in `R/repertoire/sections.R`. Pipeline 01 adds
differential clonal abundance against a within-subject baseline (which a
longitudinal bulk design has and a cross-sectional panel usually does not);
pipeline 02 adds the immunoglobulin steps when the locus is `IGH`.

---

## Configuration

Copy `config/config.example.yml` to `config/config.yml` and edit it.
`config/config.yml` is git-ignored, so cohort labels and local paths stay
out of the repository history.

The template is commented throughout: each parameter says what it does and,
where the choice matters, what the consequence of changing it is. Validation
runs before any data is read, and an unexpected sample or condition label is
an error rather than a silent `NA` group — a mistyped label caught at startup
costs a second, while the same mistake caught by a downstream join survives
into a figure.

Relative paths resolve against the repository root, so a pipeline behaves the
same however it was launched.

---

## The pipelines

### 01 — Bulk repertoire (TCR or BCR)

AIRR Rearrangement TSV per sample, from MiXCR or IgBLAST. Set `locus` to
`TRB` for bulk TCR or `IGH` for bulk BCR; the unit of observation is the
molecule (UMI).

Composition QC; clonotype-definition sensitivity across the four definitions
in common use; repertoire structure and clonal space homeostasis; diversity
at native and rarefied depth; V/J gene usage at clonotype and molecule level;
overlap with each subject's own baseline; differential clonal abundance on
both axes (descriptive category and Fisher/FDR significance, kept separate)
with the abundance-floor sweep; convergence against a permutation null;
paired comparisons across conditions with the attainable-p bound.

This is the pipeline exercised by the example dataset.

### 02 — Single-cell repertoire (TCR or BCR)

Cell Ranger VDJ `filtered_contig_annotations.csv`, or AIRR TSV per library
(MakeDb output). The unit of observation is the cell, and the pipeline
refuses input that carries no per-cell barcode rather than counting each
molecule as a cell.

It runs the same shared sections as pipeline 01 — structure, diversity at
common cell number, gene usage, overlap, convergence — reported in cells.
With `locus: IGH` it additionally runs the immunoglobulin steps: one heavy
chain per cell with doublets removed rather than resolved; clonal assignment
within subject, threshold estimated or fixed; germlines with the D segment
masked; light-chain subclones; somatic hypermutation by IMGT region;
isotypes with the unswitched compartment reported combined. With a TCR
locus those steps are skipped, since TCRs do not hypermutate.

### 03 — Single-cell transcriptome

Cell Ranger filtered feature-barcode matrices (`.h5`) per library.

Per-sample QC with per-sample thresholds; normalisation before merging;
a dimensionality and resolution sweep before the configured choice; the
cluster-composition check that decides whether cluster proportions can be
compared at all; cluster markers; cell cycle and project-supplied gene
signatures, with signature coverage reported so that a score computed on a
third of its genes is not read as a negative result.

---

## Installation

The bulk pipeline needs only CRAN packages:

```r
install.packages(c("yaml", "dplyr", "tidyr", "ggplot2", "scales", "stringdist"))
```

The BCR pipeline additionally needs the Immcantation stack, and the
transcriptome pipeline needs Seurat:

```r
install.packages(c("alakazam", "shazam", "scoper", "dowser"))   # pipeline 02
install.packages(c("Seurat", "SeuratObject"))                   # pipeline 04
install.packages("harmony")                                     # only if gex$integration is enabled
```

All four of alakazam, shazam, scoper and dowser are needed for pipeline 02:
`readIMGT()`, `createGermlines()` and `resolveLightChains()` live in dowser,
not in the other three. Pipeline 01 needs no repertoire package beyond the
CRAN set above.

Each pipeline checks its dependencies at startup and names the ones that are
missing rather than failing partway through. Pipelines 02 and 04 also expect
reference data — IMGT germlines, specificity databases, gene signatures —
described in [data/README.md](data/README.md).

Developed against R 4.5 with Seurat 5. The bulk pipeline is also tested on
R 4.3.

---

## Citation

If these workflows contribute to published work, please cite the repository
through `CITATION.cff`, and cite the underlying methods — Change-O, SCOPer,
SHazaM, dowser, MiXCR, scRepertoire, Seurat and the rest — by their own
papers. [docs/references.md](docs/references.md) lists them with DOIs and
PMIDs, all verified against PubMed, so that a Methods section written from
this repository points at what each tool actually does.

---

## Scope

These are analysis workflows, not a package, and not a substitute for the
upstream tools. Cell Ranger, IgBLAST, MakeDb and MiXCR run outside R; the
pipelines begin from their AIRR output.

The workflows flag rather than decide. Dominant clones, subject-dominated
clusters and low-coverage signatures are reported with the information
needed to judge them, and the judgement — which should be stated in the
resulting paper — stays with the analyst.

---

## Licence

MIT. See [LICENSE](LICENSE).
