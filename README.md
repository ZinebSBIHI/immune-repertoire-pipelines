# immunerepertoire

**Shared R package behind three pipelines for adaptive immune receptor
repertoires and single-cell gene expression.**

The analyses are organised on one axis, the kind of data, with one
repository per pipeline and this package as their common core:

| Repository | Input | Receptor | What it adds |
|---|---|---|---|
| [bulk-repertoire-pipeline](https://github.com/ZinebSBIHI/bulk-repertoire-pipeline) | one AIRR file per sample (MiXCR, IgBLAST) | TCR or BCR | depth-controlled diversity, change against each subject's own baseline, clone tracking, pharmacodynamic focusing, specificity annotation, hypermutation by isotype |
| [single-cell-repertoire-pipeline](https://github.com/ZinebSBIHI/single-cell-repertoire-pipeline) | Cell Ranger V(D)J contigs or AIRR with barcodes | TCR or BCR | one heavy chain per cell, clones, selection and lineage trees, clones shared between compartments, per-cell specificity, receptor embedding joined to the transcriptome |
| [single-cell-gex-pipeline](https://github.com/ZinebSBIHI/single-cell-gex-pipeline) | Cell Ranger count matrices | — | per-sample QC, recorded clustering sweep, confounding check, markers, contrasts between conditions, module scores |

Every analysis step is a function of this package. A pipeline only fixes the
order of the steps and passes results from one to the next, so the methods
can be read, tested and cited once, and the three pipelines cannot drift
apart in what they compute.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("ZinebSBIHI/immune-repertoire-pipelines")
```

The core depends only on CRAN packages (dplyr, ggplot2, tidyr, scales,
stringdist, yaml). Each family of methods asks for its own packages when it
runs, with a message naming the step that needs them:

| Methods | Packages |
|---|---|
| B-cell clones, germlines, hypermutation, selection, lineage | alakazam, shazam, scoper, dowser (CRAN); ggtree (Bioconductor) to draw trees |
| Single-cell transcriptome | Seurat, SeuratObject; harmony if batch correction is chosen |
| Clone tracking and V-J diagrams | ggalluvial, circlize, patchwork |
| CDR3 networks and motifs | igraph, ggseqlogo; Biostrings and DECIPHER for aligned logos |
| Receptor embedding | Trex (`remotes::install_github("BorchLab/Trex")`) |

Each pipeline repository has an `install.R` that installs what it needs.

---

## What is in the package

Functions are grouped by file prefix.

| Prefix | Contents |
|---|---|
| `core-` | configuration (one validated YAML file), inputs and sample sheets, run directories, figure theme and colour-blind-safe palettes |
| `rep-` | clonotype definitions, diversity and rarefaction to a common depth, Hill spectra, CDR3 length and physicochemistry, CDR3 networks and motif logos, clone tracking and V-J chord diagrams, the repertoire sections shared by the bulk and single-cell pipelines |
| `ana-` | differential clonal abundance against a within-subject baseline, the pharmacodynamic view of focusing and new clonotypes, convergence with a permutation null, paired statistics at small sample size, specificity annotation against public references |
| `bcr-` | one heavy chain per cell, clonal inference, germlines, somatic hypermutation, isotypes, replacement/silent ratios, BASELINe selection, lineage trees |
| `sc-` | joining a repertoire to a Seurat object, per-cell specificity crossed with clonal expansion, receptor embedding and weighted nearest neighbours |
| `gex-` | per-sample QC, reduction, clustering sweep, cluster composition, markers, marker panels, condition contrasts, cell cycle, module scores |

Every `section_*()` function writes one part of a run: its tables and
figures to the run directory, its numbers and their interpretation to a
plain-text report, and its result back to the pipeline.

---

## Design

The decisions that make repertoire work go quietly wrong are made explicit
and diagnosable rather than left as defaults. Each is argued, with the code
that implements it, in **[docs/design_notes.md](docs/design_notes.md)**:

1. sequencing depth is equalised (within subject) before richness or
   clonality is compared;
2. the abundance floor is swept and its cost reported, because it moves de
   novo calls in both directions at once;
3. descriptive category and statistical significance are separate axes, and
   no output asserts that a clonotype was induced;
4. convergence between subjects is read against a label-permutation null;
5. paired tests report the smallest p-value the sample size can attain, and
   direction concordance across subjects;
6. a database match is a hypothesis, under one stated matching rule per data
   type;
7. a cluster dominated by one subject is named before any proportion is
   compared, and batch correction is a decision, not a default.

Every run writes its resolved configuration, `sessionInfo()`, figure
legends, a table manifest and a numeric report beside its outputs, in a
timestamped directory that no later run overwrites.

---

## Tests

```r
devtools::test()      # or R CMD check on the built package
```

The tests pin each analytical decision above to a known answer, run every
section end to end on simulated data, and check that the per-cell
specificity annotation reproduces `Trex::annotateDB()` when Trex is
installed.

---

## Citation

Cite the software with [CITATION.cff](CITATION.cff), and the methods each
analysis relies on from [docs/references.md](docs/references.md), where every
entry was checked against PubMed. The central ones:

- Vander Heiden JA et al. AIRR Community standardized representations for
  annotated immune repertoires. *Front Immunol.* 2018;9:2206.
  [doi:10.3389/fimmu.2018.02206](https://doi.org/10.3389/fimmu.2018.02206)
- Gupta NT et al. Change-O: a toolkit for analyzing large-scale B cell
  immunoglobulin repertoire sequencing data. *Bioinformatics.*
  2015;31(20):3356-3358.
  [doi:10.1093/bioinformatics/btv359](https://doi.org/10.1093/bioinformatics/btv359)
- Hao Y et al. Dictionary learning for integrative, multimodal and scalable
  single-cell analysis. *Nat Biotechnol.* 2024;42(2):293-304.
  [doi:10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y)
- Shugay M et al. VDJdb: a curated database of T-cell receptor sequences
  with known antigen specificity. *Nucleic Acids Res.* 2018;46(D1):D419-D427.
  [doi:10.1093/nar/gkx760](https://doi.org/10.1093/nar/gkx760)

## License

MIT. See [LICENSE.md](LICENSE.md).
