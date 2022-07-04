#!/bin/sh

Dsuite="/home/irojas/bin/Dsuite/Build/Dsuite"
vcf="/home/irojas/cocci/data/VDC.vcf.gz"
meta="/home/irojas/cocci/meta"
out="/home/irojas/cocci/out"


$Dsuite Dtrios $vcf $meta/Dsuite_scenario_1.1.txt -t -o $out/Dsuite.scenario_1.1
$Dsuite Dtrios $vcf $meta/vcf_scenario1.2_subset.txt -o $out/Dsuite.scenario_1.2
$Dsuite Dtrios $vcf $meta/vcf_scenario1.3_subset.txt -o $out/Dsuite.scenario_1.3
$Dsuite Dtrios $vcf $meta/vcf_scenario1.4_subset.txt -o $out/Dsuite.scenario_1.4
$Dsuite Dtrios $vcf $meta/vcf_scenario1.5_subset.txt -o $out/Dsuite.scenario_1.5
$Dsuite Dtrios $vcf $meta/vcf_scenario1.6_subset.txt -o $out/Dsuite.scenario_1.6
$Dsuite Dtrios $vcf $meta/vcf_scenario1.7_subset.txt -o $out/Dsuite.scenario_1.7
$Dsuite Dtrios $vcf $meta/vcf_scenario1.8_subset.txt -o $out/Dsuite.scenario_1.8
$Dsuite Dtrios $vcf $meta/vcf_scenario1.9_subset.txt -o $out/Dsuite.scenario_1.9
$Dsuite Dtrios $vcf $meta/vcf_scenario1.10_subset.txt -t $meta/cocci_1.10.tree -o $out/Dsuite.scenario_1.10
$Dsuite Dtrios $vcf $meta/vcf_scenario1_subset.txt -o $out/Dsuite.scenario_1
$Dsuite Dtrios $vcf $meta/vcf_scenario2_subset.txt -t $meta/cocci_2.0.tree -o $out/Dsuite.scenario_2.0
$Dsuite Dtrios $vcf $meta/vcf_scenario3_subset.txt -t $meta/cocci_3.0.tree -o $out/Dsuite.scenario_3.0

