We should improve the logging from prolfquapp. 

most of it is nice. 

[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] using : /opt/r-libs-site/prolfqua
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] using : /opt/r-libs-site/prolfquapp
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] YAML file read: config.yaml
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] prolfquapp yaml


Butn then we again and again we have:
[2026-07-06 09:25:57] copy /opt/r-libs-site/prolfquapp/doc/Grp2Analysis_V2_R6.Rmd to /work/Grp2Analysis_V2_R6.Rmd
[2026-07-06 09:25:57] copy /opt/r-libs-site/prolfquapp/doc/DiffExpQC_R6.Rmd to /work/DiffExpQC_R6.Rmd
[2026-07-06 09:25:57] your working directory now should contain: 2 new files :

Stuff which is messaged instead of being logged.




```
[2026-07-06 09:25:42] → prolfqua_dea: docker run --rm --init --user 10263:666 --mount type=bind,source=/scratch/A414_DEA/WU348193/work/input,target=/work -w /work docker.io/prolfqua/prolfquapp:2.3.3 prolfqua_dea.sh -y config.yaml --software prolfquapp.DIANN -m firth --outdir C38944WU348193 --indir .
[2026-07-06 09:25:43] ++ Rscript --vanilla -e 'cat(system.file(package = '''prolfquapp'''))'
[2026-07-06 09:25:43] + PACKAGE_PATH=/opt/r-libs-site/prolfquapp
[2026-07-06 09:25:43] + R_SCRIPT_PATH=/opt/r-libs-site/prolfquapp/application/CMD_DEA_V2.R
[2026-07-06 09:25:43] + [[ -f /opt/r-libs-site/prolfquapp/application/CMD_DEA_V2.R ]]
[2026-07-06 09:25:43] + echo 'Rscript --vanilla "/opt/r-libs-site/prolfquapp/application/CMD_DEA_V2.R" "-y' config.yaml --software prolfquapp.DIANN -m firth --outdir C38944WU348193 --indir '."'
[2026-07-06 09:25:43] + Rscript --vanilla /opt/r-libs-site/prolfquapp/application/CMD_DEA_V2.R -y config.yaml --software prolfquapp.DIANN -m firth --outdir C38944WU348193 --indir .
[2026-07-06 09:25:43] Rscript --vanilla "/opt/r-libs-site/prolfquapp/application/CMD_DEA_V2.R" "-y config.yaml --software prolfquapp.DIANN -m firth --outdir C38944WU348193 --indir ."
[2026-07-06 09:25:44] Warning message:
[2026-07-06 09:25:44] package ‘optparse’ was built under R version 4.6.1
[2026-07-06 09:25:56] INFO [2026-07-06 07:25:56] LIBRARY PATHS (.libPaths()):/opt/r-libs-site
[2026-07-06 09:25:57] /usr/local/lib/R/site-library
[2026-07-06 09:25:57] /usr/lib/R/site-library
[2026-07-06 09:25:57] /usr/lib/R/library
[2026-07-06 09:25:57] Warning message:
[2026-07-06 09:25:57] package ‘prolfquapp’ was built under R version 4.6.1
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] using : /opt/r-libs-site/prolfqua
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] using : /opt/r-libs-site/prolfquapp
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] YAML file read: config.yaml
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] prolfquapp yaml
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] <list>
[2026-07-06 09:25:57] ├─indir: "."
[2026-07-06 09:25:57] ├─dataset: "dataset.csv"
[2026-07-06 09:25:57] ├─yaml: "config.yaml"
[2026-07-06 09:25:57] ├─software: "prolfquapp.DIANN"
[2026-07-06 09:25:57] ├─outdir: "C38944WU348193"
[2026-07-06 09:25:57] ├─model: "firth"
[2026-07-06 09:25:57] ├─flat_outdir: FALSE
[2026-07-06 09:25:57] └─help: FALSE
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] Writing to output directory : C38944WU348193/DEA_20260706_PI38944_O38944_WU348193_vsn and file :prolfqua_202607060725.log
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] prolfquapp paramters :
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] <list>
[2026-07-06 09:25:57] ├─ext_reader: <list>
[2026-07-06 09:25:57] │ ├─dataset<chr [0]>: ""
[2026-07-06 09:25:57] │ ├─extra_args: "list()"
[2026-07-06 09:25:57] │ ├─preprocess<chr [0]>: ""
[2026-07-06 09:25:57] │ └─get_files<chr [0]>: ""
[2026-07-06 09:25:57] ├─group: "G_"
[2026-07-06 09:25:57] ├─path: "C38944WU348193"
[2026-07-06 09:25:57] ├─flat_outdir: FALSE
[2026-07-06 09:25:57] ├─zipdir_name: "DEA_20260706_PI38944_O38944_WU34..."
[2026-07-06 09:25:57] ├─prefix: "DEA"
[2026-07-06 09:25:57] ├─software: "prolfquapp.DIANN"
[2026-07-06 09:25:57] ├─project_spec: <list>
[2026-07-06 09:25:57] │ ├─input_URL: "https://fgcz-bfabric.uzh.ch/bfab..."
[2026-07-06 09:25:57] │ ├─workunit_Id: 348193
[2026-07-06 09:25:57] │ ├─order_Id: 38944
[2026-07-06 09:25:57] │ ├─project_name: ""
[2026-07-06 09:25:57] │ └─project_Id: 38944
[2026-07-06 09:25:57] └─processing_options: <list>
[2026-07-06 09:25:57] ├─internal<chr [0]>: ""
[2026-07-06 09:25:57] ├─model: "firth"
[2026-07-06 09:25:57] ├─model_missing: TRUE
[2026-07-06 09:25:57] ├─interaction: FALSE
[2026-07-06 09:25:57] ├─nr_peptides: 1
[2026-07-06 09:25:57] ├─remove_decoys: TRUE
[2026-07-06 09:25:57] ├─remove_cont: TRUE
[2026-07-06 09:25:57] ├─FDR_threshold: 0.1
[2026-07-06 09:25:57] ├─diff_threshold: 1
[2026-07-06 09:25:57] ├─aggregate: "medpolish"
[2026-07-06 09:25:57] └─transform: "vsn"
[2026-07-06 09:25:57] copy /opt/r-libs-site/prolfquapp/doc/Grp2Analysis_V2_R6.Rmd to /work/Grp2Analysis_V2_R6.Rmd
[2026-07-06 09:25:57] copy /opt/r-libs-site/prolfquapp/doc/DiffExpQC_R6.Rmd to /work/DiffExpQC_R6.Rmd
[2026-07-06 09:25:57] your working directory now should contain: 2 new files :
[2026-07-06 09:25:57] 
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] Software: prolfquapp.DIANN
[2026-07-06 09:25:57] Rows: 58 Columns: 4
[2026-07-06 09:25:57] ── Column specification ────────────────────────────────────────────────────────
[2026-07-06 09:25:57] Delimiter: ","
[2026-07-06 09:25:57] chr (4): Relative Path, Name, Grouping Var, CONTROL
[2026-07-06 09:25:57] 
[2026-07-06 09:25:57] ℹ Use `spec()` to retrieve the full column specification for this data.
[2026-07-06 09:25:57] ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] levels: c("HSA_spleen", "HSA_Tumor") c("C", "T")
[2026-07-06 09:25:57] HSA_Tumor HSA_spleen
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] Files data: ./WU329547_report.tsv
[2026-07-06 09:25:57] INFO [2026-07-06 07:25:57] Files fasta: ./WU329547_report-fasta-database.fasta
[2026-07-06 09:27:48] Rows: 9039170 Columns: 60
[2026-07-06 09:27:48] ── Column specification ────────────────────────────────────────────────────────
[2026-07-06 09:27:48] Delimiter: "t"
[2026-07-06 09:27:48] chr (12): File.Name, Run, Protein.Group, Protein.Ids, Protein.Names, Genes, ...
[2026-07-06 09:27:48] dbl (48): PG.Quantity, PG.Normalised, PG.MaxLFQ, Genes.Quantity, Genes.Norma...
[2026-07-06 09:27:48] 
[2026-07-06 09:27:48] ℹ Use `spec()` to retrieve the full column specification for this data.
[2026-07-06 09:27:48] ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
[2026-07-06 09:48:50] INFO [2026-07-06 07:48:50] nr : 58 files annotated out of 217
[2026-07-06 09:48:51] Joining with `by = join_by(raw.file)`
[2026-07-06 09:48:53] column sampleName already exists, using :Name
[2026-07-06 09:49:39] completing cases
[2026-07-06 09:49:48] completing cases done
[2026-07-06 09:49:48] setup done
[2026-07-06 09:49:49] completing cases
[2026-07-06 09:50:01] INFO [2026-07-06 07:50:01] start reading fasta.
[2026-07-06 09:50:01] INFO [2026-07-06 07:50:01] get_annot : ./WU329547_report-fasta-database.fasta
[2026-07-06 09:50:03] INFO [2026-07-06 07:50:03] get_annot : finished reading
[2026-07-06 09:50:04] INFO [2026-07-06 07:50:04] get_annot : extract headers
[2026-07-06 09:50:04] INFO [2026-07-06 07:50:04] get_annot : all seq : 27126
[2026-07-06 09:50:04] INFO [2026-07-06 07:50:04] get_annot : isUniprot : TRUE
[2026-07-06 09:50:06] INFO [2026-07-06 07:50:06] get_annot : extracted gene names
[2026-07-06 09:50:06] INFO [2026-07-06 07:50:06] get_annot : protein length
[2026-07-06 09:50:09] INFO [2026-07-06 07:50:09] get_annot : nr of tryptic peptides per protein computed.
[2026-07-06 09:50:09] INFO [2026-07-06 07:50:09] reading fasta done, creating protein annotation.
[2026-07-06 09:50:09] INFO [2026-07-06 07:50:09] protein annotation done.
[2026-07-06 09:50:09] INFO [2026-07-06 07:50:09] start fitting / contrasts for model: firth
[2026-07-06 09:50:09] INFO [2026-07-06 07:50:09] AGGREGATING PEPTIDE DATA: medpolish.
[2026-07-06 09:50:09] Column added : log_Peptide.Quantity
[2026-07-06 09:50:11] starting aggregation
[2026-07-06 09:57:12] completing cases
[2026-07-06 09:57:13] Column added : exp_medpolish
[2026-07-06 09:57:13] INFO [2026-07-06 07:57:13] END OF PROTEIN AGGREGATION
[2026-07-06 09:57:13] INFO [2026-07-06 07:57:13] Transforming using vsn::justvsn
[2026-07-06 09:57:28] Joining with `by = join_by(Name, isotopeLabel, protein_Id)`
[2026-07-06 09:57:28] INFO [2026-07-06 07:57:28] Transforming data : vsn.
[2026-07-06 09:57:28] INFO [2026-07-06 07:57:28] model formula: normalized_abundance ~ G_
[2026-07-06 09:57:28] completing cases
[2026-07-06 09:57:29] Joining with `by = join_by(protein_Id)`
```