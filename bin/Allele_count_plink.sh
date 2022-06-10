##Count alleles for Site Frequency Spectrum
#Note the individual file, which contain the IDs of the samples to include. Two columns must be contained per sample

#Count alleles including - non-genic regions
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-OV --freq --counts --out ../data/Allele_counts/Nongenic/Cult-OV
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-SMOCC --freq --counts --out ../data/Allele_counts/Nongenic/Cult-SMOCC
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-SMOCC-BlaTla --freq --counts --out ../data/Allele_counts/Nongenic/Cult-SMOCC-BlaTla
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-SUR --freq --counts --out ../data/Allele_counts/Nongenic/Cult-SUR
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-TMVB --freq --counts --out ../data/Allele_counts/Nongenic/Cult-TMVB
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Cult-TMVB-Esp --freq --counts --out ../data/Allele_counts/Nongenic/Cult-TMVB-Esp
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Feral --freq --counts --out ../data/Allele_counts/Nongenic/Feral
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-SMOCC-EspDia --freq --counts --out ../data/Allele_counts/Nongenic/Wild-SMOCC-EspDia
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-SMOCC-Rego --freq --counts --out ../data/Allele_counts/Nongenic/Wild-SMOCC-Rego
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-striatus --freq --counts --out ../data/Allele_counts/Nongenic/Wild-striatus
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-TMVB-CdMex --freq --counts --out ../data/Allele_counts/Nongenic/Wild-TMVB-CdMex
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-TMVB-SanJoa --freq --counts --out ../data/Allele_counts/Nongenic/Wild-TMVB-SanJoa
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-TMVB-Tepoz --freq --counts --out ../data/Allele_counts/Nongenic/Wild-TMVB-Tepoz
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-SUR-CH --freq --counts --out ../data/Allele_counts/Nongenic/Wild-SUR-CH
plink --file ../data/SNP_categ/cocci_no_genic --keep ../meta/Pops/Wild-SUR-O --freq --counts --out ../data/Allele_counts/Nongenic/Wild-SUR-O

#Count alleles - introns
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-OV --freq --counts --out ../data/Allele_counts/Introns/Cult-OV
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-SMOCC --freq --counts --out ../data/Allele_counts/Introns/Cult-SMOCC
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-SMOCC-BlaTla --freq --counts --out ../data/Allele_counts/Introns/Cult-SMOCC-BlaTla
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-SUR --freq --counts --out ../data/Allele_counts/Introns/Cult-SUR
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-TMVB --freq --counts --out ../data/Allele_counts/Introns/Cult-TMVB
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Cult-TMVB-Esp --freq --counts --out ../data/Allele_counts/Introns/Cult-TMVB-Esp
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Feral --freq --counts --out ../data/Allele_counts/Introns/Feral
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-SMOCC-EspDia --freq --counts --out ../data/Allele_counts/Introns/Wild-SMOCC-EspDia
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-SMOCC-Rego --freq --counts --out ../data/Allele_counts/Introns/Wild-SMOCC-Rego
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-striatus --freq --counts --out ../data/Allele_counts/Introns/Wild-striatus
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-TMVB-CdMex --freq --counts --out ../data/Allele_counts/Introns/Wild-TMVB-CdMex
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-TMVB-SanJoa --freq --counts --out ../data/Allele_counts/Introns/Wild-TMVB-SanJoa
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-TMVB-Tepoz --freq --counts --out ../data/Allele_counts/Introns/Wild-TMVB-Tepoz
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-SUR-CH --freq --counts --out ../data/Allele_counts/Introns/Wild-SUR-CH
plink --file ../data/SNP_categ/cocci_intron --keep ../meta/Pops/Wild-SUR-O --freq --counts --out ../data/Allele_counts/Introns/Wild-SUR-O

#Count alleles including only synosymous mutations	
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-OV --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-OV-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-SMOCC --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SMOCC-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-SMOCC-BlaTla --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SMOCC-BlaTla-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-SUR --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SUR-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-TMVB --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-TMVB-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Cult-TMVB-Esp --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-TMVB-Esp-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Feral --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Feral-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-SMOCC-EspDia --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SMOCC-EspDia-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-SMOCC-Rego --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SMOCC-Rego-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-striatus --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-striatus-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-TMVB-CdMex --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-CdMex-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-TMVB-SanJoa --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-SanJoa-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-TMVB-Tepoz --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-Tepoz-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-SUR-CH --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SUR-CH-syn
plink --file ../data/coccineus --extract ../meta/SNPs_synonym --keep ../meta/Pops/Wild-SUR-O --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SUR-O-syn

#Count alleles including only nonsynosymous mutations
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-OV --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-OV-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-SMOCC --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SMOCC-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-SMOCC-BlaTla --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SMOCC-BlaTla-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-SUR --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-SUR-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-TMVB --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-TMVB-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Cult-TMVB-Esp --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Cult-TMVB-Esp-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Feral --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Feral-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-SMOCC-EspDia --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SMOCC-EspDia-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-SMOCC-Rego --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SMOCC-Rego-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-striatus --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-striatus-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-TMVB-CdMex --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-CdMex-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-TMVB-SanJoa --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-SanJoa-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-TMVB-Tepoz --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-TMVB-Tepoz-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-SUR-CH --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SUR-CH-nonsyn
plink --file ../data/coccineus --extract ../meta/SNPs_nonsynonym --keep ../meta/Pops/Wild-SUR-O --freq --counts --out ../data/Allele_counts/Syn-Nosyn/Wild-SUR-O-nonsyn
