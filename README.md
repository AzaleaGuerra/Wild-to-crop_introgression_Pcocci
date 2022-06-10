# Wild-to-crop_introgression_Pcocci

## The genomic signature of wild‐to‐crop introgression during the domestication of scarlet runner bean (*Phaseolus coccineus* L.)

This repository contains the script used in the article "The genomic signature of wild‐to‐crop introgression during the domestication of scarlet runner bean (*Phaseolus coccineus* L.)"

## Prerequistes
PLINK 1.07

VCFtools 0.1.15

Admixture 1.3

FastTreeMP

ADZE 1.0

fastsimcoal2

easySFS

### R packages
rehh

vcfR

Hierfstat

VariantAnnotation

geosphere

ggpubr

stringr

HDplot

adegenet

SNPRelate


## Directories
**bin**: contains the scripts and is the working directory for them.

**meta**: information related to data collection, populations, etc.

**data**: filtered vcf files are available at [https://osf.io/h7sa5/](https://osf.io/h7sa5/)

**fastsimcoal**: contains the script to perform the demographic analysis 

##Script organization
**1 Filtering**
1.1 Depth and missingness filtering: ``Filtering&transform.sh``
1.1 HW filtering: ``HW-filtering.sh``
1.3 Paralog detection: ``HDplot_cocci.R``
1.4 Relatedness: ``relatedness_filtering.R``

**2 Population structure**
2.1 PCA: ``PCA_coccineus.R`` and ``PCA_cocci_cult&wild.R``
2.2 Ancestry clustering: ``Admixture_coccineus.sh`` and ``Plot_admixture.R``
2.3 Phylogeny: ``FastTree_VDC.sh``

**3 Genetic diversity and differentiation**
3.1 Het, FIS and Fst: ``He-Fis-Fst_Hierfstat.R``
3.2 Privite alleles identification: ``private_alleles_pop.R`` and ``private_alleles_status.R``
3.3 Distance and He correlation: ``correlation_dist-He.R``

**4 Gene flow and introgression**


**5 Selection**
5.1 Selective sweeps: ``rehh_selection_test.R``

**6 Demography**
6.1 Site Frequency Sprectrum: ``Allele_count_plink.sh``
6.2 Runs of Homozygosity: ``ROH.sh`` and ``Plot_ROH.R``
6.3 Demographic scenarios: markdown file in the directory ``fastsimcoal`` contains further details about the construction and test of the demographic scenarios. 

## Contact
Azalea Guerrra Garcia

azalea.guerra@iecologia.unam.mx

azalea.guerra.g@gmail.com
