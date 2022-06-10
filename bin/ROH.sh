
##Estimate long runs of homoztgosity 
plink --file ../data/coccineus --homozyg-kb 300 --out ../data/LRHomoz/cocci_300kb

plink --file ../data/coccineus --homozyg-kb 500 --out ../data/LRHomoz/cocci_500kb
	
plink --file ../data/coccineus --homozyg-kb 1000 --out ../data/LRHomoz/cocci_1000kb
	