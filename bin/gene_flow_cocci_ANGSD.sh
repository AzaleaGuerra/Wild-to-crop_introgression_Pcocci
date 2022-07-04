#!/bin/bash

home="~/cocci/"
bin="../bin"
meta="../meta"
out="../out"
angsd="/home/irojas/bin/angsd/"
#Gene flow models
# Note: ANGSD  uses as input a list that contains the path to the bam files

#Models
for i in 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.1 3.1 4.1 4.1 ; do #

$angsd/angsd -doAbbababa2 1 -bam $meta/cocci_${i}.subset_bamlist -sizeFile $meta/cocci_${i}.subset_popsize -doCounts 1  -#out $out/cocci_${i}_subset -useLast 1 -nThreads 3 -minQ 20 -minMapQ 30;

Rscript $angsd/R/estAvgError.R angsdFile="$out/cocci_${i}_subset" sizeFile=$meta/cocci_$i.subset_popsize nameFile=$meta/cocci_$i.popnames > $out/Model_cocci_$i.subset.Dst.TerminalOut;
done


