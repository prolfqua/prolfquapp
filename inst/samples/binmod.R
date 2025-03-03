library(prolfqua)
library(tidyverse)
pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 5,with_missing = TRUE,weight_missing = 0.4)
dd <-LFQData$new(pep$data,pep$config)
dd$complete_cases()
dd$hierarchy_counts()

xd <- dd$data
xd
xd <- xd |> dplyr::mutate(binresp = ifelse(is.na(abundance), 0, 1))
x <- xd |> group_by(protein_Id) |> nest()
xd <- x$data[[1]] |> select(group_, binresp, peptide_Id )

head(xd)

f <- formula(binresp ~ group_ + peptide_Id)
tt <- ftable(f, xd)
multiplier = 1
offset = 0.5
tt <- tt * multiplier + offset
DFT <- as.data.frame(tt)
head(DFT)
DFT$binresp
DFT$Freq
glm( f ,
     data = DFT ,
     weights = Freq,
     family = stats::binomial)






#############
install.packages("logistf")
library(logistf)
f <- formula(binresp ~ group_ + peptide_Id)

# Aggregate the data using ftable
tt <- ftable(f, xd)

# Set multiplier and offset to adjust counts (to avoid perfect separation)

# Fit the logistic regression model using Firth's penalized likelihood

logistf(f, data = DFT, weights = Freq)

# View the summary of the fitted model
summary(modelFirth)


