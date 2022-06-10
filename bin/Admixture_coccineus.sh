#! /bin/bash

###All samples of P. coccineus
cd ../admixture_linux-1.3.0/

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; \
do ./admixture --cv ../data/coccineus.bed $K | tee ../data/Admixture/log${K}.out;
done

grep -h CV ../data/Admixture/log*.out > ../data/Admixture/coccineus_Kerror.txt

#Subset P. coccineus
cd ../admixture_linux-1.3.0/

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do ./admixture --cv ../data/cocci_subset.bed $K | tee ../data/Admixture_subset/log${K}.out;
done

grep -h CV ../data/Admixture_subset/log*.out > ../data/Admixture_subset/coccineus_Kerror.txt
