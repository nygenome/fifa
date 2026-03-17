#! /usr/bin/Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.r-project.org")
}
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(mobster))

args = commandArgs(trailingOnly=TRUE)
FFPE <- args[1]
vcf_path <- args[2]
path_out <- args[3]

sample <- FFPE
vcf <- VariantAnnotation::readVcf(vcf_path)

vcf_data <- data.frame(
    chrom = as.character(seqnames(rowRanges(vcf))),
    pos = start(rowRanges(vcf)),
    REF = ref(vcf),
    ALT = alt(vcf),  
    VAF = unlist(geno(vcf)$AF[,sample], use.names=FALSE))%>%
dplyr::select(-c(ALT.group, ALT.group_name)) %>%
dplyr::rename(ALT = ALT.value) %>%
as_tibble()
        
fit <- vcf_data %>%
    dplyr::filter(VAF >= 0.05 & VAF < 1) %>%
    mobster_fit(
    .,   
    K = c(1,2,3), 
    samples = 1, 
    init = 'random',
    tail = TRUE,
    epsilon = 1e-06,
    maxIter = 100, 
    fit.type = 'MM', 
    seed = 12345,
    model.selection = 'reICL',
    trace = FALSE,
    parallel = FALSE,
    pi_cutoff = 0.02,
    N_cutoff = 10
)

sample_data <- vcf_data %>%
    dplyr::left_join(Clusters(fit$best), by=c('chrom', 'pos', 'REF', 'ALT')) %>% 
    dplyr::select(c(chrom,pos,REF,ALT,Tail)) %>%
    dplyr::mutate(Tail = tidyr::replace_na(Tail, 1),
                    sample=sample) %>%
    dplyr::select(c(sample,chrom,pos,REF,ALT,Tail))

write.csv(sample_data, file=path_out, row.names=FALSE)

