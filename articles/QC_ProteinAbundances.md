# Protein Signal Intensities within Groups

- Order :
- Workunit :

![Number of proteins quantified per
sample.](QC_ProteinAbundances_files/figure-html/hierachycountspersample-1.png)

Number of proteins quantified per sample.

The columns of the table below contain:

- protein_Id - protein identifier
- nrPeptides - number of peptides assigned to the protein
- description - protein description from the FASTA header
- nrMeasured\_ - number of samples in which the protein was quantified,
  by group and overall
- meanAbundance\_ - mean protein abundance per group
- signal_percent\_ - percentage of the total protein signal attributed
  to the protein

(ref:proteinCumulative) Protein-level signal contribution. Each point
represents one protein. Proteins are ordered on the x-axis by their
contribution to the total signal; the y-axis shows the percentage of the
total iBAQ (intensity-based absolute quantification) signal attributed
to each protein.

The iBAQ signal per protein is calculated by summing the
`Precursor.Quantity` values across all precursors assigned to the
protein and dividing by the number of theoretically observable peptides.
A large share of the recorded signal can be concentrated in a small set
of highly abundant proteins. If the dominant proteins are the cleavage
enzyme, common contaminants such as human keratins, or the bait protein,
this can indicate shortcomings in sample preparation or cleanup.

(ref:proteinCumulative)

(ref:nrProtplots) Mean protein abundance per group (column
meanAbundance\_) as a function of the number of peptides (column
nrPeptides).
