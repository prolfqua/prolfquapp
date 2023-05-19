library(SummarizedExperiment)
xx <- readRDS("deResult.rds")
xx
names(metadata(xx))
names(metadata(xx)$param)
assays(xx)

rownames(assay(xx))

colnames(rowData(xx))


# type - protein
# biotypes - "protein"
# gc - NA
# GO BP - NA
# GO CC - NA
# entriz_id - if possible
# featWidth - protein length
# isPresentProble = usedInTest => TRUE FALSE model Type.
