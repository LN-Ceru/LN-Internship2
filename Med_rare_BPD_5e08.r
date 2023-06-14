#Ella loads in proteins, and ID linker file and polygenic scores

######merging the datasets#########
setwd("/exports/igmm/eddie/GenScotDepression/users/Laurence/assoc_input/lmem_msters/")
norm_prot_read = read.csv("./GS_mediumrare_InvNorm.txt", header = TRUE, sep = '')
PRS_5e08 =  read.csv("./BPD_PRS_5e08.score", header = TRUE, sep = '')
covars = read.csv("./GS_covar_mediumrare_10.txt", header = TRUE, sep = '')

#merging on the PRS for MDD first

norm_prot_prs = merge(norm_prot_read, PRS_5e08, by.x = 'IID', by.y = 'IID')
norm_prot_id = merge(norm_prot_prs, covars[,1:14], by.x = "IID", by.y = "IID")

#extracting the marker names for the protein columns
markers <- norm_prot_read[,c(3:4237)]
marker_names <- colnames(markers)

### LOAD IN COXME REQUIREMENTS
.libPaths("/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(survival,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(quadprog,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(kinship2,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(coxme,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(readxl,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
#library(tidyverse,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
library(gdata,lib.loc = "/exports/igmm/eddie/GenScotDepression/users/Laurence/software")
#load in the pedigree file with the family data in

ped <- read.csv("./pedigree_formatted.csv")

#create a kinship matrix for GS

# Create kinship matrix for GS
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid)) # Pedigree list with 26 total subjects in 5 families
kin_model <- kinship(kin)

# Function to Extract Lmekin Results
extract_coxme_table <- function (mod){
  #beta <- mod$coefficients #$fixed is not needed
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- beta/se
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

#make a vector length == the number of proteins
length <- 4235

#initialise the results table with the columns the size for all proteins (=length)
norm_prot_results <- data.frame(SeqId = 1:length, bpdPRS_beta = 1:length, bpdPRS_SE = 1:length, bpdPRS_P = 1:length)



#single protein test
#prot_name <- as.character(colnames(markers[,1]))
#norm_mod <- lmekin(scale(norm_prot_id[,prot_name]) ~ scale(norm_prot_id$Pt_5e.08) + as.factor(norm_prot_id$sex) + scale(norm_prot_id$age) + (1|norm_prot_id$IID), data = norm_prot_id, varlist = kin_model*2)

#print(prot_name)
#norm_prot_results[,1] = prot_name
#norm_prot_results[,2] <- extract_coxme_table(norm_mod)[2,1]
#norm_prot_results[,3] <- extract_coxme_table(norm_mod)[2,2]
#norm_prot_results[,4] <- extract_coxme_table(norm_mod)[2,4]

#for loop for running the lmekin model for all proteins
for(i in 1:4235){
  prot_name <- as.character(colnames(markers[i]))

  norm_mod <- lmekin(scale(norm_prot_id[,prot_name]) ~ scale(norm_prot_id$Pt_5e.08) + as.factor(norm_prot_id$sex) + scale(norm_prot_id$age) + scale(norm_prot_id$C1) + scale(norm_prot_id$C2) + scale(norm_prot_id$C3) + scale(norm_prot_id$C4) + scale(norm_prot_id$C5) + scale(norm_prot_id$C6) + scale(norm_prot_id$C7) + scale(norm_prot_id$C8) + scale(norm_prot_id$C9) + scale(norm_prot_id$C10) + (1|norm_prot_id$IID), data = norm_prot_id, varlist = kin_model*2)

  print(i)
  print(prot_name)
  norm_prot_results[i,1] <- prot_name
  norm_prot_results[i,2] <- extract_coxme_table(norm_mod)[2,1]
  norm_prot_results[i,3] <- extract_coxme_table(norm_mod)[2,2]
  norm_prot_results[i,4] <- extract_coxme_table(norm_mod)[2,4]

}

# save a copy of results tables change to EDDIE locations:
write.csv(norm_prot_results, "/exports/igmm/eddie/GenScotDepression/users/Laurence/assoc_output/BPD_mediumrare_prot_res.csv", row.names = F, quote = F)
