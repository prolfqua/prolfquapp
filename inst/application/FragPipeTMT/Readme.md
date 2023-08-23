# FAQ

## How the protein abundances are estimated

We use the information stored in the psm.tsv to derive the protein abundance estimates. In brief, we filter the psm's to have a
Purity > 0.5, and a PeptideProphetProb > 0.9.
Afterwards, we add all the assigned psm's for each modified peptide sequence. Finally, to estimate the protein abundances, we use the Tukeys median polish method, which uses the peptide intensities to give robust estimate of the protein abundance in each sample. (The protein abundance is approximately (since Tukeys median polish does some additional estimation) the median abundance of all the peptides of a protein in the sample.

## How is the number of peptides per portein deterimined.

The numer of peptides per protein information is taken from the psm.tsv file (has a column nrPeptides).

## Where can I find the psm.tsv file

The psm.tsv file is the primary output of the FragPipe workflow. More information about each of the columns in the psm.tsv file is provided here:
https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html#psmtsv


