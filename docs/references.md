# References

Every reference below was checked against PubMed; DOIs and PMIDs are
reproduced as that record gives them. Entries without a PMID are not indexed
in PubMed and are cited from the source named.

Software is cited by its method paper rather than by its package name, so
that a Methods section written from this repository points at what the tool
actually does.

## Standards and formats

- Vander Heiden JA, Marquez S, Marthandan N, et al. AIRR Community
  standardized representations for annotated immune repertoires.
  *Front Immunol.* 2018;9:2206.
  [doi:10.3389/fimmu.2018.02206](https://doi.org/10.3389/fimmu.2018.02206).
  PMID 30323809.
  — The Rearrangement schema every pipeline here reads and writes.

- Lefranc MP, Pommié C, Ruiz M, et al. IMGT unique numbering for
  immunoglobulin and T cell receptor variable domains and Ig superfamily
  V-like domains. *Dev Comp Immunol.* 2003;27(1):55-77.
  [doi:10.1016/s0145-305x(02)00039-3](https://doi.org/10.1016/s0145-305x(02)00039-3).
  PMID 12477501.
  — The region definitions behind `IMGT_V` and `IMGT_VDJ`.

## Repertoire assembly and error correction

- Bolotin DA, Poslavsky S, Mitrophanov I, et al. MiXCR: software for
  comprehensive adaptive immunity profiling. *Nat Methods.*
  2015;12(5):380-381.
  [doi:10.1038/nmeth.3364](https://doi.org/10.1038/nmeth.3364).
  PMID 25924071.

- Ye J, Ma N, Madden TL, Ostell JM. IgBLAST: an immunoglobulin variable
  domain sequence analysis tool. *Nucleic Acids Res.* 2013;41(Web Server
  issue):W34-W40.
  [doi:10.1093/nar/gkt382](https://doi.org/10.1093/nar/gkt382).
  PMID 23671333.

- Shugay M, Britanova OV, Merzlyak EM, et al. Towards error-free profiling
  of immune repertoires. *Nat Methods.* 2014;11(6):653-655.
  [doi:10.1038/nmeth.2960](https://doi.org/10.1038/nmeth.2960).
  PMID 24793455.
  — Why `repertoire$min_abundance` is counted in UMIs and not in reads.

- Egorov ES, Merzlyak EM, Shelenkov AA, et al. Quantitative profiling of
  immune repertoires for minor lymphocyte counts using unique molecular
  identifiers. *J Immunol.* 2015;194(12):6155-6163.
  [doi:10.4049/jimmunol.1500215](https://doi.org/10.4049/jimmunol.1500215).
  PMID 25957172.

## Clonal inference and somatic hypermutation

- Gupta NT, Vander Heiden JA, Uduman M, et al. Change-O: a toolkit for
  analyzing large-scale B cell immunoglobulin repertoire sequencing data.
  *Bioinformatics.* 2015;31(20):3356-3358.
  [doi:10.1093/bioinformatics/btv359](https://doi.org/10.1093/bioinformatics/btv359).
  PMID 26069265.

- Gupta NT, Adams KD, Briggs AW, Timberlake SC, Vigneault F, Kleinstein SH.
  Hierarchical clustering can identify B cell clones with high confidence in
  Ig repertoire sequencing data. *J Immunol.* 2017;198(6):2489-2499.
  [doi:10.4049/jimmunol.1601850](https://doi.org/10.4049/jimmunol.1601850).
  PMID 28179494.
  — The source of the single-linkage default in `bcr$clone_linkage`. The
  paper reports specificity, sensitivity and positive predictive value above
  99% for single linkage, with a significant loss of sensitivity under other
  linkages, and finds that distance metrics incorporating hypermutation bias
  do not outperform plain Hamming distance.

- Nouri N, Kleinstein SH. Optimized threshold inference for partitioning of
  clones from high-throughput B cell repertoire sequencing data.
  *Front Immunol.* 2018;9:1687.
  [doi:10.3389/fimmu.2018.01687](https://doi.org/10.3389/fimmu.2018.01687).
  PMID 30093903.
  — The mixture-model threshold estimate used when
  `bcr$clone_threshold_mode: estimate`.

- Yaari G, Vander Heiden JA, Uduman M, et al. Models of somatic hypermutation
  targeting and substitution based on synonymous mutations from
  high-throughput immunoglobulin sequencing data. *Front Immunol.*
  2013;4:358.
  [doi:10.3389/fimmu.2013.00358](https://doi.org/10.3389/fimmu.2013.00358).
  PMID 24298272.

- Yaari G, Uduman M, Kleinstein SH. Quantifying selection in high-throughput
  immunoglobulin sequencing data sets. *Nucleic Acids Res.* 2012;40(17):e134.
  [doi:10.1093/nar/gks457](https://doi.org/10.1093/nar/gks457).
  PMID 22641856.
  — BASELINe.

- Hoehn KB, Pybus OG, Kleinstein SH. Phylogenetic analysis of migration,
  differentiation, and class switching in B cells. *PLoS Comput Biol.*
  2022;18(4):e1009885.
  [doi:10.1371/journal.pcbi.1009885](https://doi.org/10.1371/journal.pcbi.1009885).
  PMID 35468128.
  — dowser: lineage trees (`section_lineage_trees()`).

- Hershberg U, Uduman M, Shlomchik MJ, Kleinstein SH. Improved methods for
  detecting selection by mutation analysis of Ig V region sequences.
  *Int Immunol.* 2008;20(5):683-694.
  [doi:10.1093/intimm/dxn026](https://doi.org/10.1093/intimm/dxn026).
  PMID 18397909.
  — Why replacement and silent mutations are counted separately in CDR and
  FWR (`section_rs_ratio()`).

- Uduman M, Yaari G, Hershberg U, Stern JA, Shlomchik MJ, Kleinstein SH.
  Detecting selection in immunoglobulin sequences. *Nucleic Acids Res.*
  2011;39(Web Server issue):W499-W504.
  [doi:10.1093/nar/gkr413](https://doi.org/10.1093/nar/gkr413).
  PMID 21665923.
  — The focused test statistic used by `section_baseline_selection()`.

- Schliep KP. phangorn: phylogenetic analysis in R. *Bioinformatics.*
  2011;27(4):592-593.
  [doi:10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706).
  PMID 21169378.
  — The maximum-likelihood tree engine behind dowser's `pml` build.

## Single-cell repertoire and transcriptome

- Borcherding N, Bormann NL, Kraus G. scRepertoire: an R-based toolkit for
  single-cell immune receptor analysis. *F1000Res.* 2020;9:47.
  [doi:10.12688/f1000research.22139.2](https://doi.org/10.12688/f1000research.22139.2).
  PMID 32789006.

- Yang Q, Safina KR, Nguyen KDQ, Tuong ZK, Borcherding N. scRepertoire 2:
  enhanced and efficient toolkit for single-cell immune profiling.
  *PLoS Comput Biol.* 2025;21(6):e1012760.
  [doi:10.1371/journal.pcbi.1012760](https://doi.org/10.1371/journal.pcbi.1012760).
  PMID 40577285.

- Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative,
  multimodal and scalable single-cell analysis. *Nat Biotechnol.*
  2024;42(2):293-304.
  [doi:10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y).
  PMID 37231261.
  — Seurat v5, including the layered objects that make `JoinLayers()`
  necessary before differential expression.

- Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal
  single-cell data. *Cell.* 2021;184(13):3573-3587.e29.
  [doi:10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048).
  PMID 34062119.
  — Seurat v4; weighted-nearest-neighbour multimodal integration.

- Stuart T, Butler A, Hoffman P, et al. Comprehensive integration of
  single-cell data. *Cell.* 2019;177(7):1888-1902.e21.
  [doi:10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031).
  PMID 31178118.
  — `vst` variable feature selection, `LogNormalize`, `ScaleData`, `RunPCA`,
  `FindNeighbors`, `FindClusters`.

- Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate
  integration of single-cell data with Harmony. *Nat Methods.*
  2019;16(12):1289-1296.
  [doi:10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0).
  PMID 31740819.
  — Used only when `gex$integration$method: harmony`, which is not the
  default.

- Becht E, McInnes L, Healy J, et al. Dimensionality reduction for
  visualizing single-cell data using UMAP. *Nat Biotechnol.* 2019;37:38-44.
  [doi:10.1038/nbt.4314](https://doi.org/10.1038/nbt.4314).
  PMID 30531897.

- Borcherding N, Kim W, Quinn M, et al. CD4+ T cells exhibit distinct
  transcriptional phenotypes in the lymph nodes and blood following mRNA
  vaccination in humans. *Nat Immunol.* 2024;25(9):1731-1741.
  [doi:10.1038/s41590-024-01888-9](https://doi.org/10.1038/s41590-024-01888-9).
  PMID 39164479.
  — Trex: receptor embedding and the bundled reference of receptors of
  known specificity, whose annotation rule `annotate_cells_by_cdr3()`
  reproduces.

- Schattgen SA, Guion K, Crawford JC, et al. Integrating T cell receptor
  sequences and transcriptional profiles by clonotype neighbor graph analysis
  (CoNGA). *Nat Biotechnol.* 2022;40(1):54-63.
  [doi:10.1038/s41587-021-00989-2](https://doi.org/10.1038/s41587-021-00989-2).
  PMID 34426704.
  — One representative cell per clonotype (`conga_reduction()`).

- Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular
  ecosystem of metastatic melanoma by single-cell RNA-seq. *Science.*
  2016;352(6282):189-196.
  [doi:10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501).
  PMID 27124452.
  — The cell cycle gene sets behind `cc.genes.updated.2019`, and the module
  score approach used by `AddModuleScore`.

- Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E. Fast unfolding of
  communities in large networks. *J Stat Mech.* 2008;2008(10):P10008.
  [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008).
  — The Louvain algorithm behind `FindClusters`. Not indexed in PubMed.

## Convergence and public clonotypes

- Soto C, Bombardi RG, Branchizio A, et al. High frequency of shared
  clonotypes in human B cell receptor repertoires. *Nature.*
  2019;566(7744):398-402.
  [doi:10.1038/s41586-019-0934-8](https://doi.org/10.1038/s41586-019-0934-8).
  PMID 30760926.
  — Essential context for any convergence claim: shared clonotypes are
  common between healthy individuals, and were found in cord blood as well
  as adult repertoires, so sharing is not by itself evidence of a shared
  antigen exposure. This is the reason `permute_public_groups()` exists.

## Differential clonal abundance

- Sahin U, Schmidt M, Derhovanessian E, et al. Individualized mRNA vaccines
  evoke durable T cell immunity in adjuvant TNBC. *Nature.*
  2026;651(8107):1088-1096.
  [doi:10.1038/s41586-025-10004-2](https://doi.org/10.1038/s41586-025-10004-2).
  PMID 41708868.
  — The source of the de novo / pre-existing distinction used by
  `test_clonal_abundance()`. **Note what is and is not taken from it.** In
  that study, candidate selection combined bulk tracking with paired
  single-cell TRA/TRB sequencing and functional validation, and only the two
  categories were used. The additional categories in this implementation
  (`stable`, `contracting`, `lost`), the fold-change threshold, the
  Fisher/FDR significance axis and the `significantly_expanded` label are
  this repository's own and must not be attributed to that paper. In
  particular, nothing here supports a "vaccine-induced" or
  "treatment-induced" label: that claim required their functional
  validation, which a bulk repertoire alone cannot provide.

- Weber JS, Carlino MS, Khattak A, et al. Individualised neoantigen therapy
  mRNA-4157 (V940) plus pembrolizumab versus pembrolizumab monotherapy in
  resected melanoma (KEYNOTE-942): a randomised, phase 2b study. *Lancet.*
  2024;403(10427):632-644.
  [doi:10.1016/S0140-6736(23)02268-7](https://doi.org/10.1016/S0140-6736(23)02268-7).
  PMID 38246194.
  — Reference for the pharmacodynamic presentation of repertoire focusing
  over time.

## Specificity databases

- Shugay M, Bagaev DV, Zvyagin IV, et al. VDJdb: a curated database of
  T-cell receptor sequences with known antigen specificity. *Nucleic Acids
  Res.* 2018;46(D1):D419-D427.
  [doi:10.1093/nar/gkx760](https://doi.org/10.1093/nar/gkx760).
  PMID 28977646.

- Bagaev DV, Vroomans RMA, Samir J, et al. VDJdb in 2019: database
  extension, new analysis infrastructure and a T-cell receptor motif
  compendium. *Nucleic Acids Res.* 2020;48(D1):D1057-D1062.
  [doi:10.1093/nar/gkz874](https://doi.org/10.1093/nar/gkz874).
  PMID 31588507.

- Tickotsky N, Sagiv T, Prilusky J, Shifrut E, Friedman N. McPAS-TCR: a
  manually curated catalogue of pathology-associated T cell receptor
  sequences. *Bioinformatics.* 2017;33(18):2924-2929.
  [doi:10.1093/bioinformatics/btx286](https://doi.org/10.1093/bioinformatics/btx286).
  PMID 28481982.

- Goncharov M, Bagaev D, Shcherbinin D, et al. VDJdb in the pandemic era: a
  compendium of T cell receptors specific for SARS-CoV-2. *Nat Methods.*
  2022;19(9):1017-1019.
  [doi:10.1038/s41592-022-01578-0](https://doi.org/10.1038/s41592-022-01578-0).
  PMID 35970936.

- Raybould MIJ, Kovaltsuk A, Marks C, Deane CM. CoV-AbDab: the coronavirus
  antibody database. *Bioinformatics.* 2021;37(5):734-735.
  [doi:10.1093/bioinformatics/btaa739](https://doi.org/10.1093/bioinformatics/btaa739).
  PMID 32805021.
  — The layout read by `read_antibody_reference()`.

- Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB):
  2018 update. *Nucleic Acids Res.* 2019;47(D1):D339-D343.
  [doi:10.1093/nar/gky1006](https://doi.org/10.1093/nar/gky1006).
  PMID 30357391.

A caution that applies to all four: a CDR3 matching a database entry is a
hypothesis about specificity, not a measurement of it. Database coverage is
heavily skewed towards a few well-studied epitopes, so the annotated
fraction of a repertoire reflects what has been studied at least as much as
what the repertoire recognises. An unannotated clonotype is uncharacterised,
not non-specific.

## Diversity spectra, CDR3 features and clone tracking

- Hill MO. Diversity and evenness: a unifying notation and its consequences.
  *Ecology.* 1973;54(2):427-432.
  — Hill numbers (`hill_numbers()`). Not indexed in PubMed; cited from the
  journal.

- Hurlbert SH. The nonconcept of species diversity: a critique and
  alternative parameters. *Ecology.* 1971;52(4):577-586.
  — Expected richness at a smaller sample size, the exact rarefaction curve
  of `rarefaction_curve()`. Not indexed in PubMed; cited from the journal.

- Greiff V, Bhat P, Cook SC, Menzel U, Kang W, Reddy ST. A bioinformatic
  framework for immune repertoire diversity profiling enables detection of
  immunological status. *Genome Med.* 2015;7(1):49.
  [doi:10.1186/s13073-015-0169-8](https://doi.org/10.1186/s13073-015-0169-8).
  PMID 26140055.
  — Diversity profiles across Hill orders rather than a single index.

- Kyte J, Doolittle RF. A simple method for displaying the hydropathic
  character of a protein. *J Mol Biol.* 1982;157(1):105-132.
  [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0).
  PMID 7108955.
  — The hydropathy scale behind the CDR3 GRAVY score.

- Wagih O. ggseqlogo: a versatile R package for drawing sequence logos.
  *Bioinformatics.* 2017;33(22):3645-3647.
  [doi:10.1093/bioinformatics/btx469](https://doi.org/10.1093/bioinformatics/btx469).
  PMID 29036507.

- Wright ES. DECIPHER: harnessing local sequence context to improve protein
  multiple sequence alignment. *BMC Bioinformatics.* 2015;16:322.
  [doi:10.1186/s12859-015-0749-z](https://doi.org/10.1186/s12859-015-0749-z).
  PMID 26445311.
  — Alignment of the CDR3s of a public group before its logo is drawn.

- Csardi G, Nepusz T. The igraph software package for complex network
  research. *InterJournal Complex Systems.* 2006:1695.
  — CDR3 similarity networks. Not indexed in PubMed.

- Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and enhances
  circular visualization in R. *Bioinformatics.* 2014;30(19):2811-2812.
  [doi:10.1093/bioinformatics/btu393](https://doi.org/10.1093/bioinformatics/btu393).
  PMID 24930139.
  — V-J chord diagrams.

- Brunson JC. ggalluvial: layered grammar for alluvial plots.
  *J Open Source Softw.* 2020;5(49):2017.
  [doi:10.21105/joss.02017](https://doi.org/10.21105/joss.02017).
  PMID 36919162.
  — Clone tracking across samples and timepoints.

## Example data

- Laserson U, Vigneault F, Gadala-Maria D, et al. High-resolution antibody
  dynamics of vaccine-induced immune responses. *Proc Natl Acad Sci U S A.*
  2014;111(13):4928-4933.
  [doi:10.1073/pnas.1323862111](https://doi.org/10.1073/pnas.1323862111).
  PMID 24639495.
  — Source of `alakazam::ExampleDb`, the public example of the bulk B-cell
  pipeline.

## Statistics and metrics

- Chao A. Nonparametric estimation of the number of classes in a population.
  *Scand J Stat.* 1984;11(4):265-270.
  — Chao1. Not indexed in PubMed; cited from the journal.

- Wilcoxon F. Individual comparisons by ranking methods. *Biometrics Bull.*
  1945;1(6):80-83.
  [doi:10.2307/3001968](https://doi.org/10.2307/3001968).
  — Not indexed in PubMed.

- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
  and powerful approach to multiple testing. *J R Stat Soc Series B.*
  1995;57(1):289-300.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
  — Not indexed in PubMed.

- van der Loo MPJ. The stringdist package for approximate string matching.
  *R J.* 2014;6(1):111-122.
  [doi:10.32614/RJ-2014-011](https://doi.org/10.32614/RJ-2014-011).
  — Levenshtein distance. Not indexed in PubMed.

## Figures

- Okabe M, Ito K. Color universal design (CUD): how to make figures and
  presentations that are friendly to colorblind people. 2008.
  <https://jfly.uni-koeln.de/color/>
  — The categorical palette in `R/core-plotting.R`. Not a journal article.
