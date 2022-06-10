####################################################################
###################### FILTERING ###################################

##Depth
#Get site mean depth
vcftools --vcf ../data/Prefiltering/Phaseolus.v2.vcf --site-mean-depth --keep ../meta/Prefiltering/IDs_cocci_NOfilt.txt --out ../data/Prefiltering/cocci_site_depth

#Filter site, min depth = 6, excluding loci not mapped in chromosomes considering VDC group samples
vcftools --vcf ../data/Prefiltering/Phaseolus.v2.vcf --min-meanDP 6 --keep ../meta/Prefiltering/IDs_vdc_NOfilt.txt --chr Chr01 --chr Chr02 --chr Chr03 --chr Chr04 --chr Chr05 --chr Chr06 --chr Chr07 --chr Chr08 --chr Chr09 --chr Chr10 --chr Chr11 --recode --out ../data/Prefiltering/vdc_depth_filt

##Missingness on individual basis
vcftools --vcf ../data/Prefiltering/vdc_depth_filt.recode.vcf --missing-indv --out ../data/Prefiltering/missing_vdc

##Exclude samples with missingness greater than 0.3 and samples from locations with few samples (Ameca, Oaxac2001, zaca, zina, cguz)
vcftools --vcf ../data/Prefiltering/vdc_depth_filt.recode.vcf --remove ../meta/Prefiltering/high_miss-exclud.txt --min-meanDP 6 --recode --out ../data/Prefiltering/vdc_samp_filt

################ Phaseolus coccineus ###############################
##Minimum allele count = 1, max missing = 0.95 (5% of missing data), keep only biallelic SNPs
vcftools --vcf ../data/Prefiltering/vdc_samp_filt.recode.vcf --keep ../meta/IDs_cocci.txt --mac 1 --max-missing 0.95 --max-alleles 2 --recode --out ../data/Prefiltering/ccocci_site_filt

##Separate samples according to their status (wild, feral, cultivated). These files will be used for paralog identification.
vcftools --vcf ../data/Prefiltering/cocci_site_filt.vcf --keep ../meta/IDs_wild.txt --recode --out ../data/Prefiltering/wild_site_filt
vcftools --vcf ../data/Prefiltering/cocci_site_filt.vcf --keep ../meta/IDs_cult.txt --recode --out ../data/Prefiltering/cult_site_filt
vcftools --vcf ../data/Prefiltering/cocci_site_filt.vcf --keep ../meta/IDs_fer.txt --recode --out ../data/Prefiltering/fer_site_filt

##Paralogs filtering. Putative paralogs were identified using HDPlot on R
vcftools --vcf ../data/Prefiltering/cocci_site_filt.vcf --exclude ../meta/SNPs_paral_all.txt --recode --out ../data/Prefiltering/cocci_paralog_filt

## HW filtering: SNPs that were not in HW (p<0.01) in at least one wild population were excluded. Identification of these SNPs was perfomred with Plink 
vcftools --vcf ../data/Prefiltering/cocci_paralog_filt.recode.vcf --exclude ../meta/snp_filt_HW.txt --recode --out ../data/coccineus

##Tranform from vcf to plink ped
vcftools --vcf ../data/coccineus.recode.vcf --plink --out ../data/coccineus
#Subset
vcftools --vcf ../data/coccineus.recode.vcf --plink --remove ../meta/Excluded_subset.txt --out ../data/cocci_subset
./plink --file ../data/cocci_subset --make-bed --out ../data/cocci_subset

#split file separating according to their status
vcftools --vcf ../data/coccineus.recode.vcf --keep ../meta/IDs_wild.txt --recode --out ../data/wild
vcftools --vcf ../data/coccineus.recode.vcf --keep ../meta/IDs_cult.txt --recode --out ../data/cult
vcftools --vcf ../data/coccineus.recode.vcf --keep ../meta/IDs_fer.txt --recode --out ../data/fer

######################### VDC GROUP #######################################
##Minimum allele count = 1, max missing = 0.95 (5% of missing data), keep only biallelic SNPs
vcftools --vcf ../data/Prefiltering/vdc_samp_filt.recode.vcf --mac 1 --max-missing 0.95 --max-alleles 2 --recode --out ../data/Prefiltering/VDC_site_filt

##To get SNP IDs and change chromosome names, transformation of the resulting files was made in Tassel. Tassel outputs: 
../data/Prefiltering/VDC_site_filt.vcf

##Paralogs and HW filtering
vcftools --vcf ../data/Prefiltering/VDC_site_filt.vcf --exclude ../meta/SNPs_paral-HW.txt --recode --out ../data/VDC

##Tranform from vcf to plink ped using vcftools (DO NOT use Tassel)
vcftools --vcf ../data/VDC.recode.vcf --plink --out ../data/VDC



