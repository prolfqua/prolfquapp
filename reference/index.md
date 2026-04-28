# Package index

## All functions

- [`AnnotationProcessor`](https://prolfqua.github.io/prolfquapp/reference/AnnotationProcessor.md)
  : AnnotationProcessor
- [`DEAReportGenerator`](https://prolfqua.github.io/prolfquapp/reference/DEAReportGenerator.md)
  : DEAReportGenerator
- [`DEAnalyse`](https://prolfqua.github.io/prolfquapp/reference/DEAnalyse.md)
  : Differential expression analysis engine using prolfqua facade
  classes
- [`ExternalReader`](https://prolfqua.github.io/prolfquapp/reference/ExternalReader.md)
  : external reader R6 class for handling external data sources
- [`FragPipe`](https://prolfqua.github.io/prolfquapp/reference/FragPipe.md)
  : Methods for reading Fragpipe outputs
- [`LFQData_from_anndata()`](https://prolfqua.github.io/prolfquapp/reference/LFQData_from_anndata.md)
  : Convert AnnData back to LFQData + ProteinAnnotation
- [`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md)
  : Methods for reading MaxQuant outputs
- [`ProcessingOptions`](https://prolfqua.github.io/prolfquapp/reference/ProcessingOptions.md)
  : processing options R6 class
- [`ProjectSpec`](https://prolfqua.github.io/prolfquapp/reference/ProjectSpec.md)
  : project specification R6 class
- [`ProjectStructure`](https://prolfqua.github.io/prolfquapp/reference/ProjectStructure.md)
  : keep track of folder paths and create them if needed
- [`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md)
  : R6 class representing ProlfquApp configuration
- [`ProteinAnnotation`](https://prolfqua.github.io/prolfquapp/reference/ProteinAnnotation.md)
  : Decorates LFQData with a row annotation and some protein specific
  functions.
- [`ProteinDataPrep`](https://prolfqua.github.io/prolfquapp/reference/ProteinDataPrep.md)
  : ProteinDataPrep
- [`QC_generator`](https://prolfqua.github.io/prolfquapp/reference/QC_generator.md)
  : QC_generator
- [`add_RevCon()`](https://prolfqua.github.io/prolfquapp/reference/add_RevCon.md)
  : add REV and zz entries - used for testing
- [`add_contrasts_vec()`](https://prolfqua.github.io/prolfquapp/reference/add_contrasts_vec.md)
  : add vector of contrasts to annotation data frame
- [`aggregate_data()`](https://prolfqua.github.io/prolfquapp/reference/aggregate_data.md)
  : dataset transform data
- [`bfabric_url_builder()`](https://prolfqua.github.io/prolfquapp/reference/bfabric_url_builder.md)
  : build bfabric urls
- [`build_protein_annot()`](https://prolfqua.github.io/prolfquapp/reference/build_protein_annot.md)
  : build Dataset protein annot, defaults are compatible with DIANN
- [`capture_output()`](https://prolfqua.github.io/prolfquapp/reference/capture_output.md)
  : capture output of function to send it to log
- [`column_to_rownames()`](https://prolfqua.github.io/prolfquapp/reference/column_to_rownames.md)
  : convert tibble to data.frame with rownames
- [`compute_IBAQ_values()`](https://prolfqua.github.io/prolfquapp/reference/compute_IBAQ_values.md)
  : compute IBAQ values
- [`copy_DEA_R6_Files()`](https://prolfqua.github.io/prolfquapp/reference/copy_DEA_R6_Files.md)
  : copy Markdown templates for DEA (R6-based CMD_DEA_V2.R)
- [`copy_docker_script()`](https://prolfqua.github.io/prolfquapp/reference/copy_docker_script.md)
  : copy dockerfile to run the DEA app
- [`copy_shell_script()`](https://prolfqua.github.io/prolfquapp/reference/copy_shell_script.md)
  : copy shellscript to run the DEA app
- [`dataset_get_functions()`](https://prolfqua.github.io/prolfquapp/reference/dataset_get_functions.md)
  : get functions for creating datasets
- [`dataset_protein_annot()`](https://prolfqua.github.io/prolfquapp/reference/dataset_protein_annot.md)
  : Dataset protein annot
- [`dataset_template_BGS()`](https://prolfqua.github.io/prolfquapp/reference/dataset_template_BGS.md)
  : create templte dataset for BGS data
- [`dataset_template_FP_TMT()`](https://prolfqua.github.io/prolfquapp/reference/dataset_template_FP_TMT.md)
  : get dataset annotation template
- [`dataset_template_MAXQUANT()`](https://prolfqua.github.io/prolfquapp/reference/dataset_template_MAXQUANT.md)
  : create template dataset for MAXQUANT data
- [`dataset_template_MSSTATS()`](https://prolfqua.github.io/prolfquapp/reference/dataset_template_MSSTATS.md)
  : create dataset template from MSStats data.
- [`dataset_template_diann()`](https://prolfqua.github.io/prolfquapp/reference/dataset_template_diann.md)
  : create dataset template from DIANN outputs
- [`diann_output_to_peptide()`](https://prolfqua.github.io/prolfquapp/reference/diann_output_to_peptide.md)
  : Create peptide level (stripped sequences) report by aggregating
  Precursor abundances.
- [`diann_read_output()`](https://prolfqua.github.io/prolfquapp/reference/diann_read_output.md)
  : read DiaNN diann-output.tsv file
- [`example_deanalyse()`](https://prolfqua.github.io/prolfquapp/reference/example_deanalyse.md)
  : Create an example DEAnalyse object from simulated data
- [`example_qc_generator()`](https://prolfqua.github.io/prolfquapp/reference/example_qc_generator.md)
  : Create an example QC_generator object from simulated data
- [`exp2()`](https://prolfqua.github.io/prolfquapp/reference/exp2.md) :
  transform lfq data with x^2 - apply if non log data is needed
- [`extract_GN()`](https://prolfqua.github.io/prolfquapp/reference/extract_GN.md)
  : extract gene names from uniprot 1sp fasta.headers
- [`extract_contrasts()`](https://prolfqua.github.io/prolfquapp/reference/extract_contrasts.md)
  : extract contrast from annotation file
- [`feature_annotation_collapse_to_single_row()`](https://prolfqua.github.io/prolfquapp/reference/feature_annotation_collapse_to_single_row.md)
  : get best feature annotation.
- [`feature_annotation_get_best_score()`](https://prolfqua.github.io/prolfquapp/reference/feature_annotation_get_best_score.md)
  : get best feature annotation.
- [`get_BGS_files()`](https://prolfqua.github.io/prolfquapp/reference/get_BGS_files.md)
  : get BGS and fasta file location in folder
- [`get_DIANN_files()`](https://prolfqua.github.io/prolfquapp/reference/get_DIANN_files.md)
  : get report.tsv and fasta file location in folder
- [`get_FP_PSM_files()`](https://prolfqua.github.io/prolfquapp/reference/get_FP_PSM_files.md)
  : get psm.tsv and fasta file location in folder
- [`get_FP_multiSite_files()`](https://prolfqua.github.io/prolfquapp/reference/get_FP_multiSite_files.md)
  : get report.tsv and fasta file location in folder
- [`get_MQ_peptide_files()`](https://prolfqua.github.io/prolfquapp/reference/get_MQ_peptide_files.md)
  : get petpide.txt and fasta file location in folder
- [`get_MSstats_files()`](https://prolfqua.github.io/prolfquapp/reference/get_MSstats_files.md)
  : get petpide.txt and fasta file location in folder
- [`get_SIM_files()`](https://prolfqua.github.io/prolfquapp/reference/get_SIM_files.md)
  : Discover simulated data files (no real files needed)
- [`get_annot_from_fasta()`](https://prolfqua.github.io/prolfquapp/reference/get_annot_from_fasta.md)
  : get_annot_from_fasta
- [`get_config()`](https://prolfqua.github.io/prolfquapp/reference/get_config.md)
  : get configuration from yaml file or create default configuration
- [`get_dummy_files()`](https://prolfqua.github.io/prolfquapp/reference/get_dummy_files.md)
  : get_dummy_files
- [`get_mzMine_files()`](https://prolfqua.github.io/prolfquapp/reference/get_mzMine_files.md)
  : get mzmine fliles
- [`get_procfuncs()`](https://prolfqua.github.io/prolfquapp/reference/get_procfuncs.md)
  : Get all processing functions from all packages
- [`list_to_R6_app_config()`](https://prolfqua.github.io/prolfquapp/reference/list_to_R6_app_config.md)
  : read minimal yaml and convert to R6 object
- [`make_DEA_config_R6()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_config_R6.md)
  : create GRP2 configuration for differential expression analysis Use
  this function if there is no Yaml Input.
- [`make_annotated_experiment()`](https://prolfqua.github.io/prolfquapp/reference/make_annotated_experiment.md)
  : make lfqdata with row annotation
- [`make_feature_annotation()`](https://prolfqua.github.io/prolfquapp/reference/make_feature_annotation.md)
  : get feature annotation.
- [`massage_CD()`](https://prolfqua.github.io/prolfquapp/reference/massage_CD.md)
  : massage CD output compound table.
- [`normalize_path()`](https://prolfqua.github.io/prolfquapp/reference/normalize_path.md)
  : Function to normalize paths for both Windows and Linux
- [`nr_tryptic_peptides()`](https://prolfqua.github.io/prolfquapp/reference/nr_tryptic_peptides.md)
  : Compute number of tryptic peptides
- [`plot_abundance_vs_percent()`](https://prolfqua.github.io/prolfquapp/reference/plot_abundance_vs_percent.md)
  : Plot relative protein abundance as a function of rank by abundance
- [`preprocess_BGS()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_BGS.md)
  : preprocess DIANN ouput, filter by q_value and nr_peptides
- [`preprocess_CD()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_CD.md)
  : load compound discoverer (CD) files
- [`preprocess_DIANN()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_DIANN.md)
  : preprocess DIANN ouput, filter by q_value and nr_peptides
- [`preprocess_DIANN_anndata()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_DIANN_anndata.md)
  : Preprocess DIANN output and return AnnData
- [`preprocess_FP_PSM()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_FP_PSM.md)
  : preprocess FP psm, filter by purity_threshold and PeptideProphetProb
- [`preprocess_MQ_peptide()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_MQ_peptide.md)
  : preprocess MQ peptide
- [`preprocess_MSstats()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_MSstats.md)
  : preprocess MSstats file coming from FragPipe
- [`preprocess_MSstats_FPDIA()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_MSstats_FPDIA.md)
  : preprocess MSstats fragpipe
- [`preprocess_SIM()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_SIM.md)
  : Preprocess simulated data
- [`preprocess_anndata_from_lfq()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_anndata_from_lfq.md)
  : Convert LFQData + ProteinAnnotation to AnnData with round-trip
  metadata
- [`preprocess_dummy()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_dummy.md)
  : preprocess_dummy
- [`preprocess_mzMine()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_mzMine.md)
  : preprocess mzMine input
- [`preprocess_software()`](https://prolfqua.github.io/prolfquapp/reference/preprocess_software.md)
  : collects preprocess methods for various software
- [`prolfqua_preprocess_functions`](https://prolfqua.github.io/prolfquapp/reference/prolfqua_preprocess_functions.md)
  : Mapping of software to string representations of functions and
  arguments Preprocess Functions Mapping
- [`read_BF_yamlR6()`](https://prolfqua.github.io/prolfquapp/reference/read_BF_yamlR6.md)
  : read yaml file and convert to R6 configuration object
- [`read_BGS()`](https://prolfqua.github.io/prolfquapp/reference/read_BGS.md)
  : get BGS and fasta file location in folder
- [`read_annotation()`](https://prolfqua.github.io/prolfquapp/reference/read_annotation.md)
  : read annotation files
- [`read_msstats()`](https://prolfqua.github.io/prolfquapp/reference/read_msstats.md)
  : read MSstats.csv files and rollup to ProteinSequence level.
- [`read_table_data()`](https://prolfqua.github.io/prolfquapp/reference/read_table_data.md)
  : read dataset file in csv, tsv or xlsx format
- [`run_contrasts_single()`](https://prolfqua.github.io/prolfquapp/reference/run_contrasts_single.md)
  : Generate contrasts for a single-factor design
- [`run_contrasts_twofactor()`](https://prolfqua.github.io/prolfquapp/reference/run_contrasts_twofactor.md)
  : Generate contrasts for a two-factor design
- [`run_dea()`](https://prolfqua.github.io/prolfquapp/reference/run_dea.md)
  : Run differential expression analysis pipeline
- [`run_make_yaml()`](https://prolfqua.github.io/prolfquapp/reference/run_make_yaml.md)
  : Generate a DEA configuration list
- [`run_qc_preprocess()`](https://prolfqua.github.io/prolfquapp/reference/run_qc_preprocess.md)
  : Preprocess quantification data for QC reporting
- [`set_lib_path()`](https://prolfqua.github.io/prolfquapp/reference/set_lib_path.md)
  : set library path with logging
- [`set_list_to_R6()`](https://prolfqua.github.io/prolfquapp/reference/set_list_to_R6.md)
  : set arguments in list config to r6obj
- [`sim_data_protAnnot()`](https://prolfqua.github.io/prolfquapp/reference/sim_data_protAnnot.md)
  : simulate peptdata and fitting protein annotation for testing
- [`strip_rownames()`](https://prolfqua.github.io/prolfquapp/reference/strip_rownames.md)
  : Strip pattern from row names of a matrix or data.frame
- [`sync_opt_config()`](https://prolfqua.github.io/prolfquapp/reference/sync_opt_config.md)
  : Synchronize opt and config
- [`tidy_FragPipe_psm()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_psm.md)
  : read psm.tsv produced by FragPipe and convert into long format
- [`tidy_FragPipe_psm_V2()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_psm_V2.md)
  : read psm.tsv produced by FragPipe and convert into long format
- [`tidy_mzMineFeatures()`](https://prolfqua.github.io/prolfquapp/reference/tidy_mzMineFeatures.md)
  : convert mzmine features to tidy table
- [`transform_lfqdata()`](https://prolfqua.github.io/prolfquapp/reference/transform_lfqdata.md)
  : Transform lfq data using robscale, vsn or log2
- [`validate_prolfquapp_anndata()`](https://prolfqua.github.io/prolfquapp/reference/validate_prolfquapp_anndata.md)
  : Validate that an AnnData has the prolfquapp uns namespace
- [`writeLinesPaired()`](https://prolfqua.github.io/prolfquapp/reference/writeLinesPaired.md)
  : nice plot for paired analysis
- [`write_annotation_file()`](https://prolfqua.github.io/prolfquapp/reference/write_annotation_file.md)
  : Write dataset to file in csv, tsv, or xlsx format
- [`write_index_html()`](https://prolfqua.github.io/prolfquapp/reference/write_index_html.md)
  : write index.html file with links to all relevant files:
- [`zipdir_name()`](https://prolfqua.github.io/prolfquapp/reference/zipdir_name.md)
  : Generate zip directory name based on project information and date
