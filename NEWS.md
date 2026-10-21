# immunerepertoire 0.2.0

The repository becomes the shared core of three pipelines, each in its own
repository (bulk-repertoire-pipeline, single-cell-repertoire-pipeline,
single-cell-gex-pipeline), and is packaged as an installable R package.

## New method families

* **Clone tracking**: dominant clonotypes compared across samples
  (`section_clone_compare()`), baseline-dominant and significantly expanded
  clonotypes followed over time (`section_clone_tracking()`), clonotypes
  shared between compartments (`section_compartment_tracking()`), V-J chord
  diagrams per sample and for the expanded clonotypes with colours fixed
  across panels (`section_vj_circos()`, `section_expanded_vj_circos()`).
* **CDR3 and pharmacodynamics**: abundance-weighted spectratype, Hill
  spectrum and exact rarefaction curves, CDR3 charge and hydrophobicity,
  nucleotide-to-amino-acid convergence, CDR3 similarity networks with a
  threshold diagnostic, motif logos of public groups, and the detection-based
  view of repertoire focusing and new clonotypes with its depth control
  (`section_pharmacodynamics()`).
* **Specificity**: readers for VDJdb (release and web export layouts),
  McPAS-TCR, the Trex reference and antibody tables such as CoV-AbDab;
  exact annotation of bulk repertoires with bystander categories
  (`section_specificity()`); per-cell annotation with propagation within the
  dataset, crossed with clonal expansion (`section_cell_specificity()`);
  public heavy chains against known antibodies (`section_antibody_lookup()`).
* **Receptor embedding**: Trex embedding, weighted nearest neighbours of RNA
  and receptor, CoNGA representatives (`section_tcr_embedding()`).
* **Advanced B-cell analyses**: hypermutation load by compartment,
  replacement/silent ratios, BASELINe selection, class-switched against
  unswitched cells, hypermutation by isotype, V gene by CDR3 length
  overviews, the largest clones, lineage trees.
* **Single-cell transcriptome** sections: per-sample QC, reduction,
  clustering sweep, cluster composition, cluster markers, marker panels,
  condition contrasts, cell cycle and module scores.

## Restored

* Abundance-floor sweep, differential abundance as a section with
  persistence of the expanded clonotypes, B-cell clonal inference for bulk and
  single-cell data, isotypes, clonal expansion thresholds, dominant clones,
  and joining a repertoire to a Seurat object.

## Changes

* Inputs: optional sample sheet (`paths$sample_sheet`), gzipped AIRR files,
  Matrix Market expression directories as well as `.h5`.
* Reports and legends are written as UTF-8 under any locale; figures have an
  explicit white background.
* Gene calls carrying an IMGT species prefix ("Homsap IGHV3-49*01 F") are
  simplified like any other.
