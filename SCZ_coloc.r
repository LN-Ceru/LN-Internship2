#SCZ_Coloc

#MDD_HPA

#testout

#MR test script 

.libPaths("/exports/igmm/eddie/GenScotDepression/users/Laurence/software")

library(TwoSampleMR)
library(dplyr)
library(magrittr)
library(readr)
library(coloc)

#BIPOLAR DISORDER#################
proteins_uniprot_list = c("Q9Y253", "P16442", "O00478", "Q29983", "Q29980")
#proteins_uniprot_list = c( "O00478", "Q29983")

#X chromosome =Q96EU7
# uniprot = Q96EU7
# start =  47199485
# end =  47204278

#LOOP SETUP 
  #results table 
  # note the following flag has been removed: protein_name =  1:length(proteins_uniprot)


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

# '/exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/Edited2_PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv'

  exposure_data = read_exposure_data(filename = '/exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/Edited2_PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv',
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "FCON",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    samplesize_col = "N",
    #min_pval = 1e-200,
    log_pval = FALSE,
    pval_col = "P",
    ncase_col = "NCAS",
    ncontrol_col = "NCON",
    chr_col = "CHR",
    pos_col = "BP")


  exposure_data$id.exposure = 'SCZ'
  exposure_data$exposure = 'SCZ'

  print('exposure read + dim:')
  print(dim(exposure_data))

  exposure_data$MAF = ifelse(exposure_data$eaf.exposure > 0.5, 1 - exposure_data$eaf.exposure, exposure_data$eaf.exposure)
  exposure_data$varbeta = (exposure_data$se.exposure * exposure_data$se.exposure)
  exposure_data$cc_ratio = max(exposure_data$ncase.exposure) / (max(exposure_data$ncase.exposure) + max(exposure_data$ncontrol.exposure))


#for each chromosome:
#subset by chrom in full datset and select clumped hits 
#for each of the top hits define a window and add the positions to a list 
#extract accompanying rsids and append to rsid list. this gives rsids for each chomosome

#/exports/eddie/scratch/v1lnisb2/pQTL_sets/subsets/CUTDOWN_SPER.txt
#paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_input/mstrs1_2SMR/", protein, "_munged.txt")
#the list is then specified when reading in the outcome data. 
# HITS ONLY NEED TO BE READ IN ONCE - ANALYSIS LOOP STARTS HERE
# Outcome data
coloc_full_output =  list()


for (i in 1:length(proteins_uniprot_list)) {
  protein <- proteins_uniprot_list[i]
  outcome_data <- read_outcome_data(
    filename = paste0("/exports/igmm/eddie/GenScotDepression/users/Laurence/mr_input/mstrs1_2SMR/", protein, "_munged.txt"),
    sep = "\t",
    snps = exposure_data$SNP,
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF1",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    samplesize_col = "N",
    pval_col = "P",
    # min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "CHR",
    pos_col = "BP"
  )

  outcome_data$outcome <- protein
  outcome_data$id.outcome <- protein

  print("outcome read + dim:")
  print(dim(outcome_data))

  # outcome data wrangling #########################

  outcome_data$MAF <- ifelse(outcome_data$eaf.outcome > 0.5, 1 - outcome_data$eaf.outcome, outcome_data$eaf.outcome)
  outcome_data$varbeta <- (outcome_data$se.outcome * outcome_data$se.outcome)


  # removing ambigous SNPs
  print('harmonising')
  #harmsdata <- harmonise_data(exposure_data, outcome_data)


  # NOTE THIS IS NOT CORRECT NEEDS REMOVE DATA ALSO
  #exposure_harms <- harmsdata[, names(harmsdata) %in% c(names(exposure_data), "remove", "varbeta.y", "MAF.y")]
  #outcome_harms <- harmsdata[, names(harmsdata) %in% c(names(outcome_data), "remove", "varbeta.x", "MAF.x")]

  exposure_harms = exposure_data[exposure_data$SNP %in% outcome_data$SNP,]
  outcome_harms = outcome_data[outcome_data$SNP %in% exposure_data$SNP,]


  # remove mismatching positions
  #exposure_harms <- exposure_harms[!harmsdata$remove, ]
  #outcome_harms <- outcome_harms[!harmsdata$remove, ]


  # SELECT TOP SNPS FROM EXPOSURE

  # filtering and clumping in the protein# exposure data
  exposure_data_topsnps = exposure_harms[exposure_harms$pval.exposure < 5e-08, ] 
  exposure_data_topsnps_maf = exposure_data_topsnps %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)

  print("exposure_data_topsnps_maf dim")
  print(dim(exposure_data_topsnps_maf))

  if (dim(exposure_data_topsnps_maf)[1] == 0) {
    exposure_data_pre_maf <- exposure_harms %>% filter(eaf.exposure < 0.99 & eaf.exposure > 0.01)
    exposure_data_topsnps_maf <- exposure_data_pre_maf[exposure_data_pre_maf$pval.exposure == min(exposure_data_pre_maf$pval.exposure), ]
    print("min p vale reset")
    print(min(exposure_data_pre_maf$pval.exposure))
  }

  exposure_data_clumped <- clump_data(exposure_data_topsnps_maf, clump_kb = 1000, clump_r2 = 0.01)

  print("clumped dim")
  print(dim(exposure_data_clumped))

  #######
  # Colocalisation datawrangling

  exposure_data_clumped$pos.exposure <- as.numeric(exposure_data_clumped$pos.exposure)

  signif_snps <- exposure_data_clumped$SNP

  rsid_list <- c()




  # INITIALISE OUTPUT

  all_chrom_output <- matrix(nrow = 0, ncol = 10)
  all_chrom_output <- as.data.frame(all_chrom_output)


  # once each protein has been read in run through all the sigificant SNPS and run coloc for each
  #for (ch in 1:22) {

    # NOTES TO SELF THE LEAD SNPS ARE NOT SHOWING UP PRESUMABLY BECUASE THEY ARE NOT ARRIVING IN THE HARMONISED DATA
    # PRESUMABLY THE POSITIONS OF THE CLUMPED DATA DONT MATCH THOSE IN THE HARMONISED DATA
    # THERE ARE 23000 ODD SNPS IN THE RAW EXPOSURE AND ONLY 599 IN THE HARMONISED DATA SO ITS PLAUSIBLE YOU ARE PULLING DIFFERENT ONES
    # CONSIDER AGAIN IN THE MORNING.
    snp_output <- matrix(nrow = 0, ncol = 10)
    snp_output <- as.data.frame(snp_output)
    names(snp_output) = c("protein", "lead snp", "CHR", "POS", "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")


    for (j in 1:length(signif_snps)) {
      print('snp number')
      print(j)
     
      #extract top SNP rowdata 
      snp_list_chr = c()
      snp1 = signif_snps[j]
      lead_snp = exposure_harms[exposure_harms$SNP == snp1, ]$SNP[1]
      lead_pos = exposure_harms[exposure_harms$SNP == snp1, ]$pos.exposure[1]
      lead_chr = exposure_harms[exposure_harms$SNP == snp1, ]$chr.exposure[1]

      #subset by chromosome 
      ch_expos <- exposure_harms[exposure_harms$chr.exposure == lead_chr, ]
      ch_output <- outcome_harms[outcome_harms$chr.outcome == lead_chr, ]

      print("chromosome")
      print(lead_chr)

      fram2 <- ch_expos %>% filter(pos.exposure >= lead_pos - 1e6 & pos.exposure <= lead_pos + 1e6)
      # snp_list_chr = unlist(append(snp_list_chr, fram2$pos.exposure))
      snp_list_chr = unique(fram2$pos.exposure)

      ch_expos_select = subset(ch_expos, ch_expos$pos.exposure %in% snp_list_chr)
      ch_output_select = subset(ch_output, ch_output$pos.outcome %in% ch_expos_select$pos.exposure)

      #harmonising for weirdo snps
      harmsdata <- harmonise_data(ch_expos_select, ch_output_select)
      ch_expos_1 <- harmsdata[, names(harmsdata) %in% c(names(exposure_data), "remove", "varbeta.y", "MAF.y")]
      ch_output_1 <- harmsdata[, names(harmsdata) %in% c(names(outcome_data), "remove", "varbeta.x", "MAF.x")]

      dataset_1 = list(snp = ch_expos_1$SNP, position = as.numeric(ch_expos_1$pos.exposure), N = as.numeric(ch_expos_1$samplesize.exposure), beta = as.numeric(ch_expos_1$beta.exposure), MAF = as.numeric(ch_expos_1$MAF.y), varbeta = as.numeric(ch_expos_1$se.exposure)^2, pvalues = as.numeric(ch_expos_1$pval.exposure), type = "cc", s = ch_expos_1$cc_ratio[1])
      dataset_2 = list(snp = ch_output_1$SNP, position = as.numeric(ch_output_1$pos.outcome), N = as.numeric(ch_output_1$samplesize.outcome), pvalues = as.numeric(ch_output_1$pval.outcome), beta = as.numeric(ch_output_1$beta.outcome), varbeta = as.numeric(ch_output_1$se.outcome)^2, MAF = as.numeric(ch_output_1$MAF.x), type = "quant", sdY = 1)
      print("number of SNPS datasets 1 and 2")
      print(length(dataset_1$snp))
      print(length(dataset_2$snp))

      print("leadsnp")
      print(lead_snp)
      print("check data 1")
      check_dataset(dataset_1)
      print("check data 2")
      check_dataset(dataset_2)

      # Analysis and storage
      coloc <- coloc.abf(dataset_1, dataset_2)
      coloc$chrom <- lead_chr

      snp_output[j, 1] <- protein
      snp_output[j, 2] <- lead_snp
      snp_output[j, 3] <- lead_chr
      snp_output[j, 4] <- lead_pos
      snp_output[j, 5] <- coloc[[1]][[1]]
      snp_output[j, 6] <- coloc[[1]][[2]]
      snp_output[j, 7] <- coloc[[1]][[3]]
      snp_output[j, 8] <- coloc[[1]][[4]]
      snp_output[j, 9] <- coloc[[1]][[5]]
      snp_output[j, 10] <- coloc[[1]][[6]]
    }
    ## WRITE RESULTS FOR PROTEIN
  
    # bind results
    coloc_full_output[[i]] =  snp_output
}

coloc_full_output <- do.call(rbind, coloc_full_output)

coloc_full_output <- as.data.frame(coloc_full_output)
coloc_full_output = coloc_full_output[rev(order(coloc_full_output$PP.H4)), ]

write.csv(coloc_full_output, "/exports/igmm/eddie/GenScotDepression/users/Laurence/coloc_output/SCZ_coloc1.csv", row.names = F)

