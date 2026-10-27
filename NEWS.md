# immunerepertoire 0.2.0

Shared R package behind three pipelines, each in its own repository
(bulk-repertoire-pipeline, single-cell-repertoire-pipeline,
single-cell-gex-pipeline). Every analysis step is a `section_*()` function
of this package; a pipeline only fixes their order.

## Repertoire

* **Structure and depth**: composition and depth QC
  (`section_composition()`), sensitivity of the clonotype count to its
  definition (`section_clonotype_definitions()`), depth, richness, clonality
  and size classes (`section_structure()`), CDR3 length and V gene usage
  (`section_gene_usage()`), pairwise overlap (`section_overlap()`).
* **Diversity**: native and common-depth diversity, rarefied within subject
  (`section_diversity()`); Hill spectrum, exact rarefaction curves and D50
  (`section_diversity_spectrum()`).
* **Change against a baseline**: abundance-floor sweep with the cost of each
  floor (`section_abundance_floor_sweep()`), differential clonal abundance
  against each subject's own baseline with persistence of the expanded
  clonotypes (`section_differential_abundance()`), convergence between
  subjects against a label-permutation null (`section_convergence()`), paired
  comparisons that report the smallest attainable p-value and direction
  concordance (`section_paired_comparisons()`).
* **Clone tracking**: dominant clonotypes compared across samples
  (`section_clone_compare()`), baseline-dominant and significantly expanded
  clonotypes followed over time (`section_clone_tracking()`), clonotypes
  shared between compartments (`section_compartment_tracking()`), V-J chord
  diagrams per sample and for the expanded clonotypes with colours fixed
  across panels (`section_vj_circos()`, `section_expanded_vj_circos()`).
* **CDR3 and pharmacodynamics**: abundance-weighted spectratype, CDR3 charge
  and hydrophobicity, nucleotide-to-amino-acid convergence
  (`section_cdr3_structure()`), features of the expanded clonotypes and
  their CDR3 similarity network with a threshold diagnostic
  (`section_cdr3_features_expanded()`), motif logos of public groups
  (`section_cdr3_logos()`), and the detection-based view of repertoire
  focusing and new clonotypes with its depth control
  (`section_pharmacodynamics()`).

## Specificity

* Readers for VDJdb (release and web export layouts), McPAS-TCR, the Trex
  reference and antibody tables such as CoV-AbDab.
* Exact annotation of bulk repertoires with bystander categories
  (`section_specificity()`) and composition of the reference entries of the
  target species (`section_reference_composition()`).
* Per-cell annotation with propagation within the dataset, crossed with
  clonal expansion (`section_cell_specificity()`).
* Public heavy chains against antibodies of known target
  (`section_antibody_lookup()`).

## B cells

* One heavy chain per cell, clonal inference for bulk and single-cell data,
  germlines and somatic hypermutation (`section_bcr_clones()`), isotypes
  (`section_isotypes()`), clonal expansion with a threshold matched to the
  material (`section_clonal_expansion()`), dominant clones
  (`section_dominant_clones()`).
* Hypermutation load by group (`section_mutation_load()`),
  replacement/silent ratios (`section_rs_ratio()`), BASELINe selection
  (`section_baseline_selection()`), class-switched against unswitched cells
  (`section_shm_switched()`), hypermutation by isotype
  (`section_mutation_by_isotype()`), IGHV by CDR3 length overview
  (`section_vgene_cdr3_overview()`), the largest clones
  (`section_top_clones()`), lineage trees (`section_lineage_trees()`).

## Single cell

* Clonotypes and clone sizes joined to a Seurat object
  (`section_seurat_attachment()`).
* Receptor embedding with Trex, weighted nearest neighbours of RNA and
  receptor, CoNGA representatives (`section_tcr_embedding()`).
* Transcriptome: per-sample QC (`section_gex_qc()`), reduction
  (`section_gex_reduction()`), clustering sweep (`section_gex_clustering()`),
  cluster composition (`section_cluster_composition()`), cluster markers
  (`section_cluster_markers()`), marker panels (`section_marker_panels()`),
  condition contrasts (`section_condition_contrasts()`), cell cycle
  (`section_cell_cycle()`) and module scores (`section_module_scores()`).

## Inputs and outputs

* Optional sample sheet (`paths$sample_sheet`), gzipped AIRR files, Matrix
  Market expression directories as well as `.h5`.
* Reports and legends are written as UTF-8 under any locale; figures have an
  explicit white background.
* Gene calls carrying an IMGT species prefix ("Homsap IGHV3-49*01 F") are
  simplified like any other.
