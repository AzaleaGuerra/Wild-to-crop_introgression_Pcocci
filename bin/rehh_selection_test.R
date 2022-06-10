library(rehh)
library(vcfR)

#Loading data
pop <- list.files(path = "../data/pops_vcf/for_rehh/", pattern = ".recode.vcf", full.names = T)
pops <- c("Cult-SMOCC", "Cult-SUR", "Cult-TMVB")

#It must be done pero chromosome or likage group

############ CHROMOSOME 1 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "1", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.1 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.1[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.1 <- lapply(ihs.1, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.1[[i]]) > 0) {
  cr.1[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 2 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "2", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.2 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.2[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.2 <- lapply(ihs.2, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.2[[i]]) > 0) {
  cr.2[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 3 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "3", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.3 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.3[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.3 <- lapply(ihs.3, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.3[[i]]) > 0) {
  cr.3[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 4
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "4", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.4 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.4[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.4 <- lapply(ihs.4, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.4[[i]]) > 0) {
    cr.4[[i]]$pop <- rep(pops[i])  
  }
}

############ CHROMOSOME 5 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "5", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.5 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.5[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.5 <- lapply(ihs.5, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.5[[i]]) > 0) {
  cr.5[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 6 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "6", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.6 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.6[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.6 <- lapply(ihs.6, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
    if (nrow(cr.6[[i]]) > 0) {
  cr.6[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 7 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "7", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.7 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.7[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.7 <- lapply(ihs.7, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.7[[i]]) > 0) {
  cr.7[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 8 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "8", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.8 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.8[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.8 <- lapply(ihs.8, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.8[[i]]) > 0) {
  cr.8[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 9
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "9", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.9 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.9[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.9 <- lapply(ihs.9, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.9[[i]]) > 0) {
  cr.9[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 10
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "10", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.10 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.10[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.10 <- lapply(ihs.10, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.10[[i]]) > 0) {
  cr.10[[i]]$pop <- rep(pops[i])
}}

############ CHROMOSOME 11 
hh <- lapply(pop, function(x)data2haplohh(hap_file = x, vcf_reader = "data.table", chr.name = "11", polarize_vcf = F))
#Genome scan - scan_hh(hh)
res.scan <- lapply(hh, function(x)scan_hh(x, polarized = F, phased = F))
#Standarized data - ihh2ihs(res.scan, freqbin = 1)
ihs.11 <- lapply(res.scan, function(x)ihh2ihs(scan = x, freqbin = 1, min_maf = 0.05, include_freq = F))
for (i in 1:length(pop)){
  ihs.11[[i]][[1]]$pop <- rep(pops[i])
}
#Identify regions with extreme scores - Candidate regions - calc_candidate_regions(ihs, threshold = 1.5.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3)
cr.11 <- lapply(ihs.11, function(x)calc_candidate_regions(x, threshold = 1.5, window_size = 25000, overlap = 12500, min_n_mrk = 3, pval = 1.3))
for (i in 1:length(pop)){
  if (nrow(cr.11[[i]]) > 0) {
  cr.11[[i]]$pop <- rep(pops[i])
}}

####### Bind all the data toghether
IHS <- numeric(0)
for (i in 1:3) {
  IHS <- rbind(IHS, do.call("rbind", c(ihs.1[[i]][1], ihs.2[[i]][1], ihs.3[[i]][1], ihs.4[[i]][1], ihs.5[[i]][1],
                                ihs.6[[i]][1], ihs.7[[i]][1], ihs.8[[i]][1], ihs.9[[i]][1], ihs.10[[i]][1], ihs.11[[i]][1])))
}
IHS <- na.omit(IHS)
write.table(IHS, file = "../data/cand_genes/iHS_results.txt", sep = "\t", quote = F, row.names = F)

thresh <- numeric(0)
for (i in pops) {
  thresh <- c(thresh, quantile(IHS[IHS$pop == i, 3], probs = 0.95, na.rm = T))
}

can.reg <- do.call("rbind", c(cr.1, cr.2, cr.3, cr.4, cr.5, cr.6, cr.7, cr.8, cr.9, cr.10, cr.11))
write.table(can.reg, file = "../data/cand_genes/candidate_regions.txt", sep = "\t", quote = F, row.names = F)

#Plots
for (i in 1:length(pop)) {
  manhattanplot(ihs.1[[i]][1], threshold =  c(-2,2), ylim = c(-2.5,2.5), pch = 20, main = pops[i])
}

plot(reg$MEAN_MRK, reg$N_MRK)
plot(reg$MEAN_EXTR_MRK, reg$N_MRK)
