#fastsimcoal2 for modelling *P. coccineus* demography

Pages 39-43 of manual for details.
Scripts mostly based on this [tutorial](https://speciationgenomics.github.io/fastsimcoal2/).

##Inputs:
- **Observed SFS:** generated with easySFS
- **Template file: **
Describes the demographic model and the parameters of interests, which is similar to a par file, but where parameters are replaced by keywords.
- **Estimation file:**
The search range of the different parameters to estimate are defined in the following est file:

**Working directory**
For steps 1, 2, 3 and 5 the wd is each model directory (``$POPS/$PREFIX``). These scripts are located within this directory.
For steps 3-1 and 4, where the distributions of the different models are compared, the wd is the pop directory (``$POPS``). These scripts are within this directory.

##Steps:
**1-** Creation of vcf subset according to the population(s) to be included in the analysis. Then run easySFS to obtain the SFS. 
``1-SFS.sh``
```
#!/bin/bash

POPS=Cult-TMVB-Spain
PREFIX=Cult-TMVB-Spain-u

vcftools --vcf ../../../data/SNP_categ/cocci_no_genic.vcf --keep ../../pops-fastsimcoal/$POPS.txt --recode --stdout | gzip -c > $POPS.vcf.gz
rm *.log

./../../easySFS/easySFS.py -i $POPS.vcf.gz -p ../../pops-fastsimcoal/$POPS.txt -a --unfolded --preview
./../../easySFS/easySFS.py -i $POPS.vcf.gz -p ../../pops-fastsimcoal/$POPS.txt -a -f --unfolded -o SFS --proj 16

cp $POPS.vcf.gz $PREFIX.vcf.gz
cp SFS/fastsimcoal2/${POPS}_DAFpop0.obs ${PREFIX}_DAFpop0.obs
```
**Note**: Run ``easySFS`` using the preview mode. Then change --proj according to the number of samples per population.

**2-** Run fastsimcol: fastsimcoal2 should not just be run once because it might not find the global optimum of the best combination of parameters. 
Run 100 times or more and select the one with the highest likelihood.
``2-fsc_runs.sh``
```
#!/bin/bash

PREFIX=Cult-TMVB-Spain-u
NCORES=2

for i in {1..100} 
do
   mkdir run$i
   cp ../../fsc26_linux64/fsc26 run$i"/"
   cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}_DAFpop0.obs run$i"/"
   cd run$i
   ./fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -n 100000 -0 -d -M -L 40 -c $NCORES
   rm fsc26
   cd ..
done

#cat run{1..5}/${PREFIX}/${PREFIX}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$8}' | sort -k 2
#Find the best run
cp ../../fsc-selectbestrun.sh ./
./fsc-selectbestrun.sh
rm fsc-selectbestrun.sh
```

**3-** Model comparison with Likelihood distributions. 
In order to infer if the models are statistically different, the likelihood distributions for each model will be constructed. This is done by running each model with the best parameter values multple times (~100 times).
The likelihoods will differ because fastsimcoal does not compute the likelihood but rather approximates it with simulations. If the ranges of likelihoods of two models overlap, it means that they do not differ significantly.

``3-Distrub_likelihoods_bests.sh``
```
#!/bin/bash

PREFIX=Cult-TMVB-Spain-u
NCORES=2

cd bestrun
cp ../../../fsc26_linux64/fsc26 ./

# create temporary obs file with name _maxL_MSFS.obs
#cp ${PREFIX}_jointMAFpop1_0.obs ${PREFIX}_maxL_jointMAFpop1_0.obs #For two populaions
cp ${PREFIX}_DAFpop0.obs ${PREFIX}_maxL_MAFpop0.obs #For one pop


# Run fastsimcoal 100 times to get the likelihood of the observed SFS under the best parameter values with 1 mio simulated SFS.
for iter in {1..100}
do
 ./fsc26 -i ${PREFIX}_maxL.par -n 100000 -0 -m -q -c $NCORES
 # Fastsimcoal will generate a new folder called ${model}_maxL and write files in there

 # collect the lhood values (Note that >> appends to the file, whereas > would overwrite it)
 sed -n '2,3p' ${PREFIX}_maxL/${PREFIX}_maxL.lhoods  >> ${PREFIX}.lhoods

 # delete the folder with results
 rm -r ${PREFIX}_maxL/
done

#Estimate AIC
./../../../calculateAIC.sh $PREFIX

rm fsc26
rm seed.txt
cd ..
```

**Notes:** change SFS files names according to the number of populations (jointMAFpop1_0.obs for two pops or DAFpop0.obs for one pop). 

**3.1-** Plot likelihoods
``3-1plot_likelihoods_best.sh``
```
R

# Read in the likelihoods
Cult_TMVB_Spain <- scan("Cult-TMVB-Spain/bestrun/Cult-TMVB-Spain.lhoods")
Cult_TMVB_Spain_u9 <- scan("Cult-TMVB-Spain-u/bestrun/Cult-TMVB-Spain-u.lhoods")

# Plot the likelihoods

jpeg(filename = "Likelihoods_scenarios.jpg", width = 780, height = 780)
par(mfrow=c(1,1))
boxplot(range = 0, Cult_TMVB_Spain, Cult_TMVB_Spain_u9, 
	ylab="Likelihood",xaxt="n")
axis(side=1,at=1:2, labels=c("Cult_TMVB_Spain","Cult_TMVB_Spain_u9"))
dev.off()
```

**4-** AIC
In order to find the best model, the likelihoods of the best run of each model should be compared. When comparing raw likelihoods, the model with more parameters will tend to result in a better fit to the data. Therefore, the [Akaike information criterium](https://en.wikipedia.org/wiki/Akaike_information_criterion) or AIC is calculated to determine if the models differ in their likelihoods accounting for the number of parameters in each model. 
Scrip ``calculateAIC.sh`` taken from [here](https://github.com/speciationgenomics/scripts).

AIC is estimated at the end of ``3-Distrub_likelihoods_bests.sh``.
To find and plot the results:
``4-AIC-sh``
```
#!/bin/bash

for i in */bestrun/*AIC
do
echo -e `basename $i`"\t"`tail -n 1 $i` >> allmodels.AIC
done
```

**5-** Bootstrap
Once the best model is known, bootstrapping is perfomed to get an interval of the parameters. 
Create subsets of the SNPs and the SFS. 
``5-Pre_bootstrap.sh``
```
#!/bin/bash

POPS=Cult-TMVB-Spain
PREFIX=Cult-TMVB-Spain
NCORES=2

cd bestrun

# Get all lines with genomic data
zgrep -v "^#" ../$PREFIX.vcf.gz > $PREFIX.allSites

# Get the header
zgrep "^#" ../$PREFIX.vcf > header

wc -l $PREFIX.allSites
#11019 Cult-TMVB-Spain.allSites

# get 100 files with 110 sites each
split -l 110 $PREFIX.allSites $PREFIX.sites.

# Generate 50 files each with randomly concatenated blocks and compute the SFS for each:
for i in {1..100}
do
  # Make a new folder for each bootstrapping iteration:
  mkdir bs$i
  cd bs$i

  # Add the header to our new bootstrapped vcf file
  cat ../header > $PREFIX.bs.$i.vcf
  # Randomly add 100 blocks
  for r in {1..100}
  do
    cat `shuf -n1 -e ../$PREFIX.sites.*` >> ${PREFIX}.bs.$i.vcf
  done
  # Compress the vcf file again
  gzip ${PREFIX}.bs.$i.vcf

  # Make an SFS from the new bootstrapped file
  # Important to change the proj according to each pop
  ./../../../../easySFS/easySFS.py -i ${PREFIX}.bs.$i.vcf.gz -p ../../../../pops-fastsimcoal/${POPS}.txt -a -f --proj 16

  # Copy the observed SFS file into this folder renaming it to match the .tpl prefix
  #cp ../${PREFIX}_jointDAFpop1_0.obs ${PREFIX}.bs.${i}_jointDAFpop1_0.obs #For two pops
  cp ../${PREFIX}_DAFpop0.obs ${PREFIX}.bs.${i}_DAFpop0.obs #For one pop
  cp ../${PREFIX}.tpl ${PREFIX}.bs.${i}.tpl
  cp ../${PREFIX}.est ${PREFIX}.bs.${i}.est

  # Say that it is finished with iteration $i
  echo bs$i" ready"

  cd ../
done

```


Then, run fsc2 under the best model 100 times with each of these bootstrapped SFS.

``5-1bootstrap_best.sh``
```
#!/bin/bash

POPS=Cult-TMVB-Spain
PREFIX=Cult-TMVB-Spain
NCORES=2

cd bestrun

#Run the parameter estimation under the best model 100 times with each of these boostrapped SFS. This would take very long.

for bs in `seq 100`
do
  cd bs$bs
  # Run fastsimcoal 100 times:
  
  for i in `seq 100`
  do
    mkdir run$i
    cd run$i
    cp ../${PREFIX}.bs.$bs* ./
    cp ../../../../../fsc26_linux64/fsc26 ./
    ./fsc26 -t ${PREFIX}.bs.$bs.tpl -e ${PREFIX}.bs.$bs.est -m -n 100000 -0 -L 40 -M -q -c $NCORES
    cd ..
  done
  # Find the best run:
  fsc-selectbestrun.sh

  cd ..
done

cd ..
```











