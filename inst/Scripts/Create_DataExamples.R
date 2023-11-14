# to run on laptop #
# if launched using RStudio
library(rstudioapi)
# find where the script is saved (to know where to save the results)
workdir = getwd()
savedir = dirname(getActiveDocumentContext()$path)
setwd(savedir)
setwd("../Data/")

## Rdata for tests
# last update, 2023/11/14 (modification to use multiple risk factors) 
library(partialLDSC)

OA_file <- system.file("data/", "OA_GEMINI.sumstats.gz", package="partialLDSC")
T2D_file <- system.file("data/", "diabetes_type_2_GEMINI.sumstats.gz", package="partialLDSC")
BPH_file <- system.file("data/", "BPH_GEMINI.sumstats.gz", package="partialLDSC")
CHD_file <- system.file("data/", "coronary_heart_GEMINI.sumstats.gz", package="partialLDSC")

BMI_file <- system.file("data/", "BMI_Yengo_2018.txt.sumstats.gz", package="partialLDSC")

# launch analysis (using default number of blocks)
A = partial_ldsc(conditions = c(OA_file, T2D_file, BPH_file, CHD_file),
                 confounders = BMI_file, 
                 condition.names = c("OA", "T2D", "BPH", "CHD"), 
                 confounder.names = "BMI",
                 ld = "~/../Documents/Exeter/Projects/Data/eur_w_ld_chr",
                 log.name = "Example_A")
saveRDS(A, file="A.RDS")


WHR_file <- system.file("data/", "whr.giant-ukbb_2018.gz.sumstats.gz", package="partialLDSC")

# launch analysis (using default number of blocks)
B = partial_ldsc(conditions = c(T2D_file, CHD_file),
                 confounders = c(WHR_file),
                 condition.names = c("T2D", "CHD"), 
                 confounder.names = c("WHR"),
                 ld = "~/../Documents/Exeter/Projects/Data/eur_w_ld_chr",
                 log.name = "Example_B")
saveRDS(B, file="B.RDS")



# launch analysis (using default number of blocks)
C = partial_ldsc(conditions = c(T2D_file, CHD_file),
                 confounders = c(BMI_file, WHR_file),
                 condition.names = c("T2D", "CHD"), 
                 confounder.names = c("BMI", "WHR"),
                 ld = "~/../Documents/Exeter/Projects/Data/eur_w_ld_chr",
                 log.name = "Example_C")

saveRDS(C, file="C.RDS")


setwd(workdir)

# make sure to "install and restart" to make sure the latest data is included in the package

