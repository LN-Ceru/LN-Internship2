#ALL_PHENOS_FORWARD_TRANS_MR
#ALL_PHENOS_FORWARD_TRANS_MR

# MARKDOWN DRAFT

# BPD_FORWARD
# BPD TRANSLOOP

# testout

# MR test script

.libPaths("/exports/igmm/eddie/GenScotDepression/users/Laurence/software")

library(TwoSampleMR)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(MungeSumstats)

# BIPOLAR DISORDER#################


proteins_uniprot_list = c("P09923", "Q9UKS6", "P49863", "O00478", "Q29983", "Q29980")
protein_start_position_list = c(232456125, 47177522, 55024256, 26440472, 31371410, 31494881)
ch_list = c(2, 11, 5, 6, 6, 6)


# LOOP SETUP
# results table
# note the following flag has been removed: protein_name =  1:length(proteins_uniprot)
RES_mr_test = data.frame(
    protein = 1:length(proteins_uniprot_list), signif_snps = 1:length(proteins_uniprot_list), min_pvalue = 1:length(proteins_uniprot_list), signif_snps_maf = 1:length(proteins_uniprot_list),
    clumped_snps = 1:length(proteins_uniprot_list), instrument = 1:length(proteins_uniprot_list), wr_effect = 1:length(proteins_uniprot_list), wr_se = 1:length(proteins_uniprot_list), wr_pval = 1:length(proteins_uniprot_list)
)


# loop for each proten
for (i in 1:length(proteins_uniprot_list)) {
    protein = proteins_uniprot_list[i]
    protein_start_pos = protein_start_position_list[i]
    ch_prot = ch_list[i]

    # mysterious bug fix - code copied and pasted straight from TwoSampleMR source code
    read_exposure_data = function(filename, clump = FALSE, sep = " ", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize", gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, chr_col = "chr", pos_col = "pos") {
        exposure_dat = data.table::fread(filename, header = TRUE, sep = sep)
        exposure_dat = format_data(
            as.data.frame(exposure_dat),
            type = "exposure",
            snps = NULL,
            phenotype_col = phenotype_col,
            snp_col = snp_col,
            beta_col = beta_col,
            se_col = se_col,
            eaf_col = eaf_col,
            effect_allele_col = effect_allele_col,
            other_allele_col = other_allele_col,
            pval_col = pval_col,
            units_col = units_col,
            ncase_col = ncase_col,
            ncontrol_col = ncontrol_col,
            samplesize_col = samplesize_col,
            gene_col = gene_col,
            id_col = id_col,
            min_pval = min_pval,
            log_pval = log_pval,
            chr_col = chr_col,
            pos_col = pos_col
        )
        exposure_dat$data_source.exposure = "textfile"
        if (clump) {
            exposure_dat = clump_data(exposure_dat)
        }
        return(exposure_dat)
    }

    # read in exposure data
    exposure_data = read_exposure_data(
        filename = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_input/mstrs1_2SMR/", protein, "_munged.txt"),
        sep = "\t",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "AF1",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        samplesize_col = "N",
        log_pval = FALSE,
        pval_col = "P",
        chr_col = "CHR",
        pos_col = "BP"
    )


    exposure_data$id.exposure = protein
    exposure_data$exposure = protein

    # read in outcome data
    outcome_data = read_outcome_data(
        filename = "/exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/pgc-bip2021-all.vcf.tsv",
        sep = "\t",
        snps = exposure_data$SNP,
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "FCON",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        pval_col = "P",
        ncase_col = "NCAS",
        ncontrol_col = "NCON",
        log_pval = FALSE,
        chr_col = "CHR",
        pos_col = "BP"
    )

    outcome_data$outcome = "BPD"
    outcome_data$id.outcome = "BPD"

    print("outcome read + dim:")
    print(dim(outcome_data))

    print("exposure read + dim:")
    print(dim(exposure_data))

    
    exposure_data = exposure_data[exposure_data$SNP %in% outcome_data$SNP,]
    outcome_data = outcome_data[outcome_data$SNP %in% exposure_data$SNP,]

    # TRANS window
    exposure_data_chr = exposure_data[exposure_data$chr.exposure == ch_prot, ]
    exposure_data_subset1 = exposure_data_chr %>% filter(pos.exposure > protein_start_pos - 1000000 & pos.exposure < protein_start_pos + 1000000)
    exposure_data_TRANS = exposure_data[!(exposure_data$SNP %in% exposure_data_subset1$SNP),]


    print("exposure_data_TRANS dim")
    print(dim(exposure_data_TRANS))


    # filtering and clumping in the protein# exposure data
    exposure_data_topsnps = exposure_data_TRANS[exposure_data_TRANS$pval.exposure < 5e-08, ]

    print("exposure_data_topsnps dim")
    print(dim(exposure_data_topsnps))

    # maf filter
    exposure_data_topsnps_maf = exposure_data_topsnps %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)


    print("exposure_data_topsnps_maf dim")
    print(dim(exposure_data_topsnps_maf))

    # if no SNPS in the TRANS window meet significance then reset the p value threshold to the lowest available
    if (dim(exposure_data_topsnps_maf)[1] == 0) {
        exposure_data_pre_maf = exposure_data_TRANS %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)
        exposure_data_topsnps_maf = exposure_data_pre_maf[exposure_data_pre_maf$pval.exposure == min(exposure_data_pre_maf$pval.exposure), ]
        print("min p vale reset")
        print(min(exposure_data_pre_maf$pval.exposure))
    }

    exposure_data_clumped = clump_data(exposure_data_topsnps_maf, clump_kb = 1000)

    # clumping and harmonising
    print("clumped dim")
    print(dim(exposure_data_clumped))


    harmonised_data = harmonise_data(exposure_dat = exposure_data_clumped, outcome_dat = outcome_data)

    print("harmonised dim dim")
    print(dim(harmonised_data))

    # using wald ratio- selecting the exposure SNP with the smallest P value
    # if there are more than 1 SNPs, then select the one with the smallest P value
    # if there is only 1 SNP, then you can just use the harmonised data directly

    if (dim(harmonised_data)[1] > 1) {
        instrumental_snpS = harmonised_data %>% filter(pval.exposure == min(pval.exposure))
        instrumental_snp = instrumental_snpS[1, ]
        mr_wr_results = mr_wald_ratio(b_exp = instrumental_snp$beta.exposure, b_out = instrumental_snp$beta.outcome, se_exp = instrumental_snp$se.exposure, se_out = instrumental_snp$se.outcome)
        print(instrumental_snp$SNP)
    } else {
        instrumental_snp = harmonised_data %>% filter(pval.exposure == min(pval.exposure))
        mr_wr_results = mr_wald_ratio(b_exp = harmonised_data$beta.exposure, b_out = harmonised_data$beta.outcome, se_exp = harmonised_data$se.exposure, se_out = harmonised_data$se.outcome)
        print(harmonised_data$SNP)
    }


    RES_mr_test[i, "protein"] = protein
    RES_mr_test[i, "signif_snps"] = dim(exposure_data_topsnps)[1]
    RES_mr_test[i, "min_pvalue"] = instrumental_snp$pval.exposure[1]
    RES_mr_test[i, "signif_snps_maf"] = dim(exposure_data_topsnps_maf)[1]
    RES_mr_test[i, "clumped_snps"] = dim(exposure_data_clumped)[1]
    RES_mr_test[i, "instrument"] = instrumental_snp$SNP[1]
    RES_mr_test[i, "wr_effect"] = mr_wr_results$b
    RES_mr_test[i, "wr_se"] = mr_wr_results$se
    RES_mr_test[i, "wr_p"] = mr_wr_results$pval

    # MR blanket test with all methods
    mr_output = mr(harmonised_data, method_list = c("mr_egger_regression", "mr_ivw", 'mr_wald_ratio'))
    assign(paste0("mr_output_", protein), mr_output)

    scatter_one <- mr_scatter_plot(mr_output, harmonised_data)
    ggsave(scatter_one[[1]], file = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_", protein, "_FORWARD_MR_scatterplot.pdf"), width = 7, height = 7)



    # get egger intercept in the event that multiple SNPs are available.
    if (dim(harmonised_data)[1] > 2) {
        egger_intercept = mr_pleiotropy_test(harmonised_data)
        write.table(egger_intercept, paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_", protein, "_FORWARD_EGGER.tsv"), sep = "\t", row.names = F, quote = F)
    }

}

# writing results for each protein #

write.table(RES_mr_test, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_prot_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)

write.table(mr_output_P09923, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_P09923_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)
write.table(mr_output_Q9UKS6, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_Q9UKS6_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)
write.table(mr_output_P49863, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_P49863_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)
write.table(mr_output_O00478, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_O00478_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)

write.table(mr_output_Q29983, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_Q29983_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)
write.table(mr_output_Q29980, "/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/BPD_TRANS_Q29980_FORWARD_WR.tsv", sep = "\t", row.names = F, quote = F)




###########################################################################################
#SCHIZOPHRENIA



rm(list=ls(all=TRUE))


###########################################################################################

#SCZ TRANSLOOP 

#testout

#MR test script 

.libPaths("/exports/igmm/eddie/GenScotDepression/users/Laurence/software")

library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(magrittr)
library(readr)
library(MungeSumstats)

#BIPOLAR DISORDER#################
proteins_uniprot_list = c("P16581","O00478")
protein_start_position_list = c(169691781, 26440700)
protein_end_position_list =   c(169703220, 26453643)
ch_list = c(1, 6)



#LOOP SETUP 
  #results table 
  # note the following flag has been removed: protein_name =  1:length(proteins_uniprot)
RES_mr_test <- data.frame(protein = 1:length(proteins_uniprot_list), signif_snps = 1:length(proteins_uniprot_list), min_pvalue = 1:length(proteins_uniprot_list), signif_snps_maf =1:length(proteins_uniprot_list) , 
clumped_snps =1:length(proteins_uniprot_list), instrument =1:length(proteins_uniprot_list), wr_effect= 1:length(proteins_uniprot_list), wr_se = 1:length(proteins_uniprot_list), wr_pval=1:length(proteins_uniprot_list))



for(i in 1:length(proteins_uniprot_list)){
  protein = proteins_uniprot_list[i]
  protein_start_pos = protein_start_position_list[i]
  protein_end_pos = protein_end_position_list[i]
  ch_prot = ch_list[i]


read_exposure_data <- function(filename, clump=FALSE, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-200, log_pval=FALSE, chr_col="chr", pos_col="pos")
{
	exposure_dat <- data.table::fread(filename, header=TRUE, sep=sep)
	exposure_dat <- format_data(
		as.data.frame(exposure_dat),
		type="exposure",
		snps=NULL,
		phenotype_col=phenotype_col,
		snp_col=snp_col,
		beta_col=beta_col,
		se_col=se_col,
		eaf_col=eaf_col,
		effect_allele_col=effect_allele_col,
		other_allele_col=other_allele_col,
		pval_col=pval_col,
		units_col=units_col,
		ncase_col=ncase_col,
		ncontrol_col=ncontrol_col,
		samplesize_col=samplesize_col,
		gene_col=gene_col,
		id_col=id_col,
		min_pval=min_pval,
		log_pval=log_pval,
		chr_col=chr_col,
		pos_col=pos_col
	)
	exposure_dat$data_source.exposure <- "textfile"
	if(clump)
	{
		exposure_dat <- clump_data(exposure_dat)
	}
	return(exposure_dat)
}


  exposure_data = read_exposure_data(filename = paste0('/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_input/mstrs1_2SMR/',protein,'_munged.txt'),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF1",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    samplesize_col = "N",
    #min_pval = 1e-200,
    log_pval = FALSE,
    pval_col = "P",
    chr_col = "CHR",
    pos_col = "BP")


  exposure_data$id.exposure = protein
  exposure_data$exposure = protein


  outcome_data = read_outcome_data(filename = '/exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/Edited2_PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv',
    sep = "\t",
    snps = exposure_data$SNP,
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "FCAS",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    ncase_col = "NCAS",
    ncontrol_col = "NCON",
    #min_pval = 1e-200,
    log_pval = FALSE,
    chr_col= "CHR",
	pos_col= "BP")

  outcome_data$outcome = 'SCZ'
  outcome_data$id.outcome = 'SCZ'

  print('outcome read + dim:')
  print(dim(outcome_data))

  print('exposure read + dim:')
  print(dim(exposure_data))

  #TRANS window 
  #exposure_data_TRANS1 = exposure_data[exposure_data$pos.exposure < protein_start_pos[i] - 500000,]  
  #exposure_data_TRANS = exposure_data[exposure_data$pos.exposure > protein_end_pos[i] + 500000,]

      
    exposure_data = exposure_data[exposure_data$SNP %in% outcome_data$SNP,]
    outcome_data = outcome_data[outcome_data$SNP %in% exposure_data$SNP,]

 # TRANS window
    exposure_data_chr = exposure_data[exposure_data$chr.exposure == ch_prot, ]
    exposure_data_subset1 = exposure_data_chr %>% filter(pos.exposure > protein_start_pos - 1000000 & pos.exposure < protein_start_pos + 1000000)
    exposure_data_TRANS = exposure_data[!(exposure_data$SNP %in% exposure_data_subset1$SNP),]


  print("exposure_data_TRANS dim")
  print(dim(exposure_data_TRANS))


  #filtering and clumping in the protein# exposure data 
  exposure_data_topsnps =  exposure_data_TRANS[exposure_data_TRANS$pval.exposure < 5e-08,]
  
  print("exposure_data_topsnps dim")
  print(dim(exposure_data_topsnps))

  exposure_data_topsnps_maf = exposure_data_topsnps %>% filter(eaf.exposure <0.99 & eaf.exposure>0.01) 

  print("exposure_data_topsnps_maf dim")
  print(dim(exposure_data_topsnps_maf))

  if (dim(exposure_data_topsnps_maf)[1] == 0) {
    exposure_data_pre_maf = exposure_data_TRANS %>% filter(eaf.exposure <0.99 & eaf.exposure>0.01) 
    exposure_data_topsnps_maf = exposure_data_pre_maf[exposure_data_pre_maf$pval.exposure == min(exposure_data_pre_maf$pval.exposure),]
    print("min p vale reset")
    print(min(exposure_data_pre_maf$pval.exposure))
  }

  exposure_data_clumped = clump_data(exposure_data_topsnps_maf, clump_kb = 1000)

  print("clumped dim")
  print(dim(exposure_data_clumped))

   
  harmonised_data = harmonise_data(exposure_dat=exposure_data_clumped, outcome_dat = outcome_data)

  print("harmonised dim dim")
  print(dim(harmonised_data))

  #using wald ratio- selecting the exposure SNP with the smallest P value 
  #if there are more than 1 SNPs, then select the one with the smallest P value
  # if there is only 1 SNP, then you can just use the harmonised data directly 

  if (dim(harmonised_data)[1] > 1) {
      instrumental_snpS <- harmonised_data %>% filter(pval.exposure == min(pval.exposure))
      instrumental_snp = instrumental_snpS[1,]
      mr_wr_results <- mr_wald_ratio(b_exp = instrumental_snp$beta.exposure, b_out = instrumental_snp$beta.outcome, se_exp = instrumental_snp$se.exposure, se_out = instrumental_snp$se.outcome)
      print(instrumental_snp$SNP)
  } else {
      instrumental_snp <- harmonised_data %>% filter(pval.exposure == min(pval.exposure))
      mr_wr_results <- mr_wald_ratio(b_exp = harmonised_data$beta.exposure, b_out = harmonised_data$beta.outcome, se_exp = harmonised_data$se.exposure, se_out = harmonised_data$se.outcome)
      print(harmonised_data$SNP)
  }

  RES_mr_test[i, 'protein'] = protein
  RES_mr_test[i, 'signif_snps'] = dim(exposure_data_topsnps)[1]
  RES_mr_test[i, 'min_pvalue'] = instrumental_snp$pval.exposure[1]
  RES_mr_test[i, 'signif_snps_maf'] = dim(exposure_data_topsnps_maf)[1]
  RES_mr_test[i, 'clumped_snps'] = dim(exposure_data_clumped)[1]
  RES_mr_test[i, 'instrument'] = instrumental_snp$SNP[1]
  RES_mr_test[i, 'wr_effect'] = mr_wr_results$b
  RES_mr_test[i, 'wr_se'] = mr_wr_results$se
  RES_mr_test[i, 'wr_p'] = mr_wr_results$pval

  mr_output = mr(harmonised_data, method_list = c("mr_egger_regression", "mr_ivw", 'mr_wald_ratio'))
  assign(paste0("mr_output_", protein), mr_output)

  scatter_one <- mr_scatter_plot(mr_output, harmonised_data)
  ggsave(scatter_one[[1]], file = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/SCZ_TRANS_", protein, "_FORWARD_MR_scatterplot.pdf"), width = 7, height = 7)


  if (dim(harmonised_data)[1] > 2) {
    egger_intercept <- mr_pleiotropy_test(harmonised_data)
    write.table(egger_intercept, paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/SCZ_", protein, "_TRANS_FORWARD_EGGER.tsv"), sep = "\t", row.names = F, quote = F)
  }
}



write.table(RES_mr_test, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/SCZ_prot_TRANS_FORWARD_WR.tsv', sep = '\t', row.names = F, quote = F)

write.table(mr_output_O00478, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/SCZ_O00478_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)
write.table(mr_output_P16581, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/SCZ_P16581_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)
 


#################################################################################
#MDD



rm(list=ls(all=TRUE))



##################################################################################



.libPaths("/exports/igmm/eddie/GenScotDepression/users/Laurence/software")

library(TwoSampleMR)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(MungeSumstats)

# BIPOLAR DISORDER#################

# extracted previously from FORWARD - automate this with a dataframe for full loop
# convert uniprot to ensembl - look into UNiprots annotation options.
# wget http://www.proteinatlas.org/<ENSEMBLE_ID_HERE>.tsv

proteins_uniprot_list = c("Q9Y253", "P16442", "O00478")
protein_start_position_list = c(43576185, 136131056, 26440472)
protein_end_position_list =   c(43620523, 136150605, 26453415)
ch_list = c(6, 9, 6)

# LOOP SETUP
# results table
# note the following flag has been removed: protein_name =  1:length(proteins_uniprot)
RES_mr_test = data.frame(
    protein = 1:length(proteins_uniprot_list), signif_snps = 1:length(proteins_uniprot_list), min_pvalue = 1:length(proteins_uniprot_list), signif_snps_maf = 1:length(proteins_uniprot_list),
    clumped_snps = 1:length(proteins_uniprot_list), instrument = 1:length(proteins_uniprot_list), wr_effect = 1:length(proteins_uniprot_list), wr_se = 1:length(proteins_uniprot_list), wr_pval = 1:length(proteins_uniprot_list)
)


# loop for each proten
for (i in 1:length(proteins_uniprot_list)) {
    protein = proteins_uniprot_list[i]
    protein_start_pos = protein_start_position_list[i]
    ch_prot = ch_list[i]

    # mysterious bug fix - code copied and pasted straight from TwoSampleMR source code
    read_exposure_data = function(filename, clump = FALSE, sep = " ", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize", gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, chr_col = "chr", pos_col = "pos") {
        exposure_dat = data.table::fread(filename, header = TRUE, sep = sep)
        exposure_dat = format_data(
            as.data.frame(exposure_dat),
            type = "exposure",
            snps = NULL,
            phenotype_col = phenotype_col,
            snp_col = snp_col,
            beta_col = beta_col,
            se_col = se_col,
            eaf_col = eaf_col,
            effect_allele_col = effect_allele_col,
            other_allele_col = other_allele_col,
            pval_col = pval_col,
            units_col = units_col,
            ncase_col = ncase_col,
            ncontrol_col = ncontrol_col,
            samplesize_col = samplesize_col,
            gene_col = gene_col,
            id_col = id_col,
            min_pval = min_pval,
            log_pval = log_pval,
            chr_col = chr_col,
            pos_col = pos_col
        )
        exposure_dat$data_source.exposure = "textfile"
        if (clump) {
            exposure_dat = clump_data(exposure_dat)
        }
        return(exposure_dat)
    }

    # read in exposure data
    exposure_data = read_exposure_data(
        filename = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_input/mstrs1_2SMR/", protein, "_munged.txt"),
        sep = "\t",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "AF1",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        samplesize_col = "N",
        log_pval = FALSE,
        pval_col = "P",
        chr_col = "CHR",
        pos_col = "BP"
    )


    exposure_data$id.exposure = protein
    exposure_data$exposure = protein

    # read in outcome data

  outcome_data = read_outcome_data(filename = '/exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/pgc-mdd2022-noGenScot-eur-v3.49.24.11.pgc',
    sep = "\t",
    snps = exposure_data$SNP,
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "FCAS",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "PVAL",
    ncase_col = "NCAS",
    ncontrol_col = "NCON",
    #min_pval = 1e-200,
    log_pval = FALSE,
    chr_col= "CHROM",
	pos_col= "POS")

    outcome_data$outcome = "MDD"
    outcome_data$id.outcome = "MDD"

    print("outcome read + dim:")
    print(dim(outcome_data))

    print("exposure read + dim:")
    print(dim(exposure_data))

    
    exposure_data = exposure_data[exposure_data$SNP %in% outcome_data$SNP,]
    outcome_data = outcome_data[outcome_data$SNP %in% exposure_data$SNP,]

    # TRANS window
    exposure_data_chr = exposure_data[exposure_data$chr.exposure == ch_prot, ]
    exposure_data_subset1 = exposure_data_chr %>% filter(pos.exposure > protein_start_pos - 1000000 & pos.exposure < protein_start_pos + 1000000)
    exposure_data_TRANS = exposure_data[!(exposure_data$SNP %in% exposure_data_subset1$SNP),]

    print("exposure_data_TRANS dim")
    print(dim(exposure_data_TRANS))


    # filtering and clumping in the protein# exposure data
    exposure_data_topsnps = exposure_data_TRANS[exposure_data_TRANS$pval.exposure < 5e-08, ]

    print("exposure_data_topsnps dim")
    print(dim(exposure_data_topsnps))

    # maf filter
    exposure_data_topsnps_maf = exposure_data_topsnps %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)


    print("exposure_data_topsnps_maf dim")
    print(dim(exposure_data_topsnps_maf))

    # if no SNPS in the TRANS window meet significance then reset the p value threshold to the lowest available
    if (dim(exposure_data_topsnps_maf)[1] == 0) {
        exposure_data_pre_maf = exposure_data_TRANS %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)
        exposure_data_topsnps_maf = exposure_data_pre_maf[exposure_data_pre_maf$pval.exposure == min(exposure_data_pre_maf$pval.exposure), ]
        print("min p vale reset")
        print(min(exposure_data_pre_maf$pval.exposure))
    }

    exposure_data_clumped = clump_data(exposure_data_topsnps_maf, clump_kb = 1000)

    # clumping and harmonising
    print("clumped dim")
    print(dim(exposure_data_clumped))


    harmonised_data = harmonise_data(exposure_dat = exposure_data_clumped, outcome_dat = outcome_data)

    print("harmonised dim dim")
    print(dim(harmonised_data))

    # using wald ratio- selecting the exposure SNP with the smallest P value
    # if there are more than 1 SNPs, then select the one with the smallest P value
    # if there is only 1 SNP, then you can just use the harmonised data directly

    if (dim(harmonised_data)[1] > 1) {
        instrumental_snpS = harmonised_data %>% filter(pval.exposure == min(pval.exposure))
        instrumental_snp = instrumental_snpS[1, ]
        mr_wr_results = mr_wald_ratio(b_exp = instrumental_snp$beta.exposure, b_out = instrumental_snp$beta.outcome, se_exp = instrumental_snp$se.exposure, se_out = instrumental_snp$se.outcome)
        print(instrumental_snp$SNP)
    } else {
        instrumental_snp = harmonised_data %>% filter(pval.exposure == min(pval.exposure))
        mr_wr_results = mr_wald_ratio(b_exp = harmonised_data$beta.exposure, b_out = harmonised_data$beta.outcome, se_exp = harmonised_data$se.exposure, se_out = harmonised_data$se.outcome)
        print(harmonised_data$SNP)
    }

    RES_mr_test[i, "protein"] = protein
    RES_mr_test[i, "signif_snps"] = dim(exposure_data_topsnps)[1]
    RES_mr_test[i, "min_pvalue"] = instrumental_snp$pval.exposure[1]
    RES_mr_test[i, "signif_snps_maf"] = dim(exposure_data_topsnps_maf)[1]
    RES_mr_test[i, "clumped_snps"] = dim(exposure_data_clumped)[1]
    RES_mr_test[i, "instrument"] = instrumental_snp$SNP[1]
    RES_mr_test[i, "wr_effect"] = mr_wr_results$b
    RES_mr_test[i, "wr_se"] = mr_wr_results$se
    RES_mr_test[i, "wr_p"] = mr_wr_results$pval

    # MR blanket test with all methods
    mr_output = mr(harmonised_data, method_list = c("mr_egger_regression", "mr_ivw", 'mr_wald_ratio'))
    assign(paste0("mr_output_", protein), mr_output)

    scatter_one <- mr_scatter_plot(mr_output, harmonised_data)
    ggsave(scatter_one[[1]], file = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_TRANS_", protein, "_FORWARD_MR_scatterplot.pdf"), width = 7, height = 7)


  if (dim(harmonised_data)[1] > 2) {
    egger_intercept <- mr_pleiotropy_test(harmonised_data)
    write.table(egger_intercept, paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_", protein, "_TRANS_FORWARD_EGGER.tsv"), sep = "\t", row.names = F, quote = F)
  }
}

write.table(RES_mr_test, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_prot_TRANS_FORWARD_WR.tsv', sep = '\t', row.names = F, quote = F)

write.table(mr_output_Q9Y253, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_Q9Y253_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)
#write.table(mr_output_Q96EU7, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_Q96EU7_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)
write.table(mr_output_P16442, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_P16442_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)
write.table(mr_output_O00478, '/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_output/cis_trans_allprots/MDD_O00478_TRANS_FORWARD_MR.tsv', sep = '\t', row.names = F, quote = F)



