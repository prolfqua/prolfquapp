# QC_generator

QC_generator

QC_generator

## Public fields

- `lfqdata`:

  lfqdata

- `lfqdata_prot`:

  lfqdata_prot

- `lfqdata_prot_IBAQ`:

  lfqdata_prot_IBAQ

- `lfqdata_prot_transformed`:

  lfqdata_prot_transformed (VSN normalized)

- `protein_annotation`:

  protein_annotation

- `lfqdata_peptide`:

  lfqdata_peptide

- `lfqdata_peptide_transformed`:

  lfqdata_peptide_transformed (VSN normalized)

- `output_dir`:

  output_dir

- `GRP2`:

  GRP2

- `TABLES2WRITE`:

  TABLES2WRITE

- `links`:

  links

## Methods

### Public methods

- [`QC_generator$new()`](#method-QC_generator-new)

- [`QC_generator$get_peptides_wide()`](#method-QC_generator-get_peptides_wide)

- [`QC_generator$get_peptides_transformed()`](#method-QC_generator-get_peptides_transformed)

- [`QC_generator$get_peptides_transformed_wide()`](#method-QC_generator-get_peptides_transformed_wide)

- [`QC_generator$get_annotation()`](#method-QC_generator-get_annotation)

- [`QC_generator$get_prot_data()`](#method-QC_generator-get_prot_data)

- [`QC_generator$get_prot_wide()`](#method-QC_generator-get_prot_wide)

- [`QC_generator$get_prot_transformed()`](#method-QC_generator-get_prot_transformed)

- [`QC_generator$get_prot_transformed_wide()`](#method-QC_generator-get_prot_transformed_wide)

- [`QC_generator$get_prot_IBAQ()`](#method-QC_generator-get_prot_IBAQ)

- [`QC_generator$get_protein_per_group_abundance()`](#method-QC_generator-get_protein_per_group_abundance)

- [`QC_generator$get_protein_per_group_abundance_with_row_annot()`](#method-QC_generator-get_protein_per_group_abundance_with_row_annot)

- [`QC_generator$get_protein_per_group_abundance_wide()`](#method-QC_generator-get_protein_per_group_abundance_wide)

- [`QC_generator$get_prot_IBAQ_wide()`](#method-QC_generator-get_prot_IBAQ_wide)

- [`QC_generator$get_list()`](#method-QC_generator-get_list)

- [`QC_generator$write_xlsx()`](#method-QC_generator-write_xlsx)

- [`QC_generator$copy_dataset()`](#method-QC_generator-copy_dataset)

- [`QC_generator$render_QC_protein_abundances()`](#method-QC_generator-render_QC_protein_abundances)

- [`QC_generator$render_sample_size_QC()`](#method-QC_generator-render_sample_size_QC)

- [`QC_generator$render_index_html()`](#method-QC_generator-render_index_html)

- [`QC_generator$render_index_md()`](#method-QC_generator-render_index_md)

- [`QC_generator$get_protein_per_group_small_wide()`](#method-QC_generator-get_protein_per_group_small_wide)

- [`QC_generator$clone()`](#method-QC_generator-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    QC_generator$new(lfqdata, protein_annotation, prolfquapp_config)

#### Arguments

- `lfqdata`:

  LFQData object

- `protein_annotation`:

  ProteinAnnotation object

- `prolfquapp_config`:

  ProlfquAppConfig object

------------------------------------------------------------------------

### Method `get_peptides_wide()`

get peptides in wide format

#### Usage

    QC_generator$get_peptides_wide()

#### Returns

peptide data in wide format

------------------------------------------------------------------------

### Method `get_peptides_transformed()`

get VSN-transformed peptide data

#### Usage

    QC_generator$get_peptides_transformed()

#### Returns

VSN-transformed peptide LFQData

------------------------------------------------------------------------

### Method `get_peptides_transformed_wide()`

get VSN-transformed peptide data in wide format

#### Usage

    QC_generator$get_peptides_transformed_wide()

#### Returns

VSN-transformed peptide data in wide format

------------------------------------------------------------------------

### Method `get_annotation()`

get annotation data

#### Usage

    QC_generator$get_annotation()

#### Returns

annotation data.frame

------------------------------------------------------------------------

### Method `get_prot_data()`

get protein data

#### Usage

    QC_generator$get_prot_data()

#### Returns

protein LFQData

------------------------------------------------------------------------

### Method `get_prot_wide()`

get protein data in wide format

#### Usage

    QC_generator$get_prot_wide()

#### Returns

protein data in wide format

------------------------------------------------------------------------

### Method `get_prot_transformed()`

get VSN-transformed protein data

#### Usage

    QC_generator$get_prot_transformed()

#### Returns

VSN-transformed protein LFQData

------------------------------------------------------------------------

### Method `get_prot_transformed_wide()`

get VSN-transformed protein data in wide format

#### Usage

    QC_generator$get_prot_transformed_wide()

#### Returns

VSN-transformed protein data in wide format

------------------------------------------------------------------------

### Method `get_prot_IBAQ()`

get IBAQ protein data

#### Usage

    QC_generator$get_prot_IBAQ()

#### Returns

IBAQ protein LFQData

------------------------------------------------------------------------

### Method `get_protein_per_group_abundance()`

get protein abundance per group

#### Usage

    QC_generator$get_protein_per_group_abundance()

#### Returns

protein abundance per group

------------------------------------------------------------------------

### Method `get_protein_per_group_abundance_with_row_annot()`

get protein abundance per group with row annotation

#### Usage

    QC_generator$get_protein_per_group_abundance_with_row_annot()

#### Returns

protein abundance per group with annotation

------------------------------------------------------------------------

### Method `get_protein_per_group_abundance_wide()`

get protein abundance per group in wide format

#### Usage

    QC_generator$get_protein_per_group_abundance_wide()

#### Returns

protein abundance per group in wide format

------------------------------------------------------------------------

### Method `get_prot_IBAQ_wide()`

get IBAQ protein data in wide format

#### Usage

    QC_generator$get_prot_IBAQ_wide()

#### Returns

IBAQ protein data in wide format

------------------------------------------------------------------------

### Method `get_list()`

get list of all tables

#### Usage

    QC_generator$get_list()

#### Returns

list of tables

------------------------------------------------------------------------

### Method `write_xlsx()`

write tables to xlsx file

#### Usage

    QC_generator$write_xlsx()

------------------------------------------------------------------------

### Method `copy_dataset()`

copy dataset/annotation file to output directory

#### Usage

    QC_generator$copy_dataset(dataset_path)

#### Arguments

- `dataset_path`:

  path to the dataset file

------------------------------------------------------------------------

### Method `render_QC_protein_abundances()`

render QC protein abundances report

#### Usage

    QC_generator$render_QC_protein_abundances()

------------------------------------------------------------------------

### Method `render_sample_size_QC()`

render sample size QC report

#### Usage

    QC_generator$render_sample_size_QC()

------------------------------------------------------------------------

### Method `render_index_html()`

render index HTML file

#### Usage

    QC_generator$render_index_html()

------------------------------------------------------------------------

### Method `render_index_md()`

render index markdown file

#### Usage

    QC_generator$render_index_md()

------------------------------------------------------------------------

### Method `get_protein_per_group_small_wide()`

get protein per group small wide format

#### Usage

    QC_generator$get_protein_per_group_small_wide()

#### Returns

protein per group data in small wide format

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    QC_generator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
