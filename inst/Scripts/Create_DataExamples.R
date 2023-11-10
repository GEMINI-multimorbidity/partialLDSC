# to run on laptop #
# if launched using RStudio
library(rstudioapi)
# find where the script is saved (to know where to save the results)
workdir = getwd()
savedir = dirname(getActiveDocumentContext()$path)
setwd(savedir)
setwd("../Data/")

## Rdata for tests
# last update, 2023/11/10 
library(partialLDSC)

OA_file <- system.file("data/", "OA_GEMINI.sumstats.gz", package="partialLDSC")
T2D_file <- system.file("data/", "diabetes_type_2_GEMINI.sumstats.gz", package="partialLDSC")
BPH_file <- system.file("data/", "BPH_GEMINI.sumstats.gz", package="partialLDSC")
BMI_file <- system.file("data/", "BMI_Yengo_2018.txt.sumstats.gz", package="partialLDSC")

# launch analysis (using default number of blocks)
A = partial_ldsc(conditions = c(OA_file, T2D_file, BPH_file),
                 confounder = BMI_file, 
                 condition.names = c("OA", "T2D", "BPH"), 
                 confounder.name = "BMI",
                 ld = "~/../Documents/Exeter/Projects/Data/eur_w_ld_chr",
                 log.name = "Example_A")
saveRDS(A, file="A.RDS")

setwd(workdir)

# make sure to "install and restart" to make sure the latest data is included in the package

