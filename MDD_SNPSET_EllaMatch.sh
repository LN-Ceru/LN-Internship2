
#!/bin/sh
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -pe sharedmem 3
#$ -l h_rt=12:00:00
#$ -M v1lnisb2@exseed.ed.ac.uk
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

Rscript ../software/PRSice.R \
    --prsice ../software/PRSice_linux \
    --base /exports/igmm/eddie/GenScotDepression/users/Laurence/sumstats_scz_bp/daner_pgc_mdd_noGenScot_eur_hg19_v3.49.24.11.neff \
    --target /exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA \
    --extract /exports/igmm/eddie/GenScotDepression/users/Laurence/Scores/MDD_ELLAMATCH_unrelated.scores.snp \
    --thread 5 \
    --stat OR \
    --pvalue P \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --lower 5e-08 \
    --maf 0.01 \
    --nonfounders \
    --base-maf FCAS:0.01,FCON:0.01 \
    --pheno /exports/igmm/eddie/GenScotDepression/users/Laurence/dum_phenos/genscot_depression_mergedFID \
    --pheno-col scid \
    --binary-target T \
    --cov-file /exports/igmm/eddie/GenScotDepression/users/Laurence/HM3mds.mds \
    --cov-col @C[1-6] \
    --no-clump \
    --bar-levels 5e-08,0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --fastscore \
    --all-score \
    --out /exports/igmm/eddie/GenScotDepression/users/Laurence/Scores/SNPSET_MDD_Ellamatch
