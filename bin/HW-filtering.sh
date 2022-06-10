
##HW filtering: p values were estimated according to populations. SNPs that were not in HW in at least one wild population were excluded. 

plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-SMOCC-EspDia --write-snplist --out ../data/HW_filt/Wild-SMOCC-EspDia
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-SMOCC-Rego --write-snplist --out ../data/HW_filt/Wild-SMOCC-Rego
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-striatus --write-snplist --out ../data/HW_filt/Wild-striatus
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-SUR --write-snplist --out ../data/HW_filt/Wild-SUR
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-TMVB-CdMex --write-snplist --out ../data/HW_filt/Wild-TMVB-CdMex
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-TMVB-SanJoa --write-snplist --out ../data/HW_filt/Wild-TMVB-SanJoa
plink --file ../data/Prefiltering/cocci_paralog_filt --hwe 0.01 --keep ../meta/Pops/Wild-TMVB-Tepoz --write-snplist --out ../data/HW_filt/Wild-TMVB-Tepoz
