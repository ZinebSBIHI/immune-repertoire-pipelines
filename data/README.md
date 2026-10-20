# Data directory

**Nothing in this directory is tracked by git, except this file and the
simulated example dataset.** `.gitignore` excludes the data formats
outright, so primary data cannot be committed by accident. Point the
pipelines at data wherever it actually lives, through `paths$data_dir` in
`config/config.yml`.

---

## Example dataset

```bash
Rscript examples/simulate_airr.R
```

writes nine simulated AIRR files to `data/example/` — three subjects, three
timepoints each — plus a `MANIFEST.csv`. These are generated from a
deliberately simple model: a log-normal baseline abundance distribution, 25
pre-existing clonotypes expanded twelve-fold at post-baseline timepoints, 20
clonotypes appearing de novo, and multinomial sampling at a depth that
varies between samples so that the depth diagnostics have something to act
on.

They carry no biological meaning and must not be used to draw any
immunological conclusion. Their purpose is to let the bulk pipeline run with
no restricted input, and to check the methods against a known answer.

---

## Expected inputs

### Pipeline 01 — single-cell TCR

Cell Ranger VDJ output, one directory per library, each containing
`filtered_contig_annotations.csv`. The library name is taken from the
containing directory.

```
data_dir/
  <library_1>/filtered_contig_annotations.csv
  <library_2>/filtered_contig_annotations.csv
```

Optionally `seurat_annotated.rds` at the top of `data_dir`, to which the
repertoire is attached. Pipeline 04 produces a suitable object; a symlink is
fine.

### Pipeline 02 — single-cell BCR

AIRR Rearrangement TSV, one file per library, as written by MakeDb.py after
IgBLAST. Searched recursively.

```
data_dir/
  <library_1>.airr.tsv
  <library_2>.airr.tsv
```

Must carry `cell_id` (or `barcode`), `locus`, `v_call`, `j_call`,
`junction`, `sequence_alignment`, `c_call`, `productive`, and a UMI count.

### Pipeline 03 — bulk repertoire

AIRR Rearrangement TSV, one file per sample, from MiXCR or IgBLAST.

Must carry `v_call`, `j_call`, `junction`, `junction_aa`, `locus`,
`productive`, and an abundance column (`umi_count` preferred,
`duplicate_count` otherwise).

### Pipeline 04 — single-cell transcriptome

Cell Ranger `count` output: `filtered_feature_bc_matrix.h5`, one per
library. Where the default file name is kept, the library name is taken from
the containing directory.

---

## File naming

Subject, timepoint and condition are parsed out of file names by the
regexes in `cohort$filename_pattern`. The pipelines stop with a readable
error listing the files that did not match, rather than producing an `NA`
group that survives into a figure.

The example dataset uses `subject_01_t0_baseline.airr.tsv`, matched by:

```yaml
filename_pattern:
  subject:   "^(subject_[0-9]+)"
  timepoint: "_(t[0-9]+)"
  condition: "_(baseline|treatment_[ab])"
```

Adjust the patterns to the existing naming rather than renaming the data;
each pattern's first capture group is what is taken, so separators written
into the pattern do not end up in the parsed value.

---

## Reference data

Set `paths$reference_dir` in the configuration. None of this is
redistributed here: each source has its own licence and its own citation
requirements, and each is versioned, so the version used belongs in Methods.

### IMGT germlines — pipeline 02, required

Fetched with Immcantation's `fetch_imgtdb.sh`, which lays out:

```
reference_dir/imgt/human/vdj/
```

read by `alakazam::readIMGT()`. Cite Lefranc et al. 2003 (PMID 12477501) and
record the download date: IMGT is updated continuously, and allele
assignments move between releases.

### Specificity databases — optional

VDJdb and McPAS-TCR, downloaded separately and placed anywhere under
`reference_dir`.

A caution worth stating in any paper that uses them: a CDR3 matching a
database entry is a hypothesis about specificity, not a measurement of one.
Coverage is heavily skewed towards a few well-studied epitopes, so the
annotated fraction of a repertoire reflects what has been studied at least
as much as what the repertoire recognises. An unannotated clonotype is
uncharacterised, not non-specific.

Cite VDJdb (PMID 28977646, 31588507) and McPAS-TCR (PMID 28481982).

### Gene signatures — pipeline 04, optional

`reference_dir/signatures.csv`, two columns:

```csv
signature,gene
plasmablast,IRF4
plasmablast,XBP1
plasmablast,MZB1
germinal_center,BCL6
germinal_center,AICDA
```

Signatures are supplied by the project rather than hard-coded, because each
one comes from a publication that should be cited where the score is
reported. The pipeline writes a coverage table giving the fraction of each
signature's genes present in the object: a signature scored on a third of
its genes is not that signature, and a low-coverage score should not be read
as a negative result.
