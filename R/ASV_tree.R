
# This script analyses dna sequance of the ASVs to create 
# a phylogenetic tree, and save it to be used in the future

#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(magrittr)
library(ape)
library(phangorn)
library(Biostrings)

rm(list=ls())

#------ run ------------

# reading the full dna sequences (fasta file)
seq_fa <- ape::read.FASTA(file="data/data_raw/data_microbiome/ASV_merged_full.fa")

# filtering for only the relevant ASVs (1951 ASVs included in the analysis)
data_asv_filtered <- read_csv("data/data_processed/data_asv_rra0.001_p0.01_th5000_all.csv") 
asv_names <- tibble(asv_ID = unique(data_asv_filtered$asv_ID))
asv_names %<>% mutate(ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", asv_ID))) %>% arrange(ID)

# filtering the fasta file and saving a new fasta file with only the relevant ASVs
seq_fa2 <- seq_fa[names(seq_fa) %in% asv_names$asv_ID]
names(seq_fa2) <- asv_names$asv_ID
ape::write.FASTA(seq_fa2, file = "data/data_raw/data_microbiome/ASV_merged_filtered.fa") # writing the filtered fasta file

# how many ASVs in the tree?
nrow(asv_names)

# align the sequences
# method 1
seq_aligned <- Biostrings::readDNAStringSet("data/data_raw/data_microbiome/ASV_merged_filtered.fa")
aligned <- DECIPHER::AlignSeqs(seq_aligned)
seq_aligned2 <- ape::as.DNAbin(aligned)

dnaphydatAll <- phangorn::phyDat(seq_aligned2, type="DNA", levels=NULL)

# finding the best tree
mt2 <- phangorn::modelTest(dnaphydatAll,
                           model="all", multicore = T, mc.cores=2)

# optimize the tree
# the function takes the best model (according to BIC) and optimizes the parameters
best_tree <- phangorn::pml_bb(mt2) 

# save the optimized tree. this is a pml object. to get the tree call: fit$tree
saveRDS(best_tree, file = "data/data_processed/phylo_tree_rra0.001_p0.01.rds")


