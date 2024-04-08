#------ 07_build_ASV_tree.r -------------------------------------
# This script analyses dna sequance of the ASVs to create 
# a phylogenetic tree, and save it to be used in the future

#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(magrittr)
library(ape)
library(adegenet)
library(phangorn)
#BiocManager::install("msa")
library(msa)

rm(list=ls())

#------ run ------------

# reading the full dna sequences (fasta file)
seq_fa <- read.FASTA(file="data/data_raw/data_microbiome/ASV_merged_full.fa")

# filtering out ASVs by relative abundance of 1%
data_asv_filtered <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv") %>% 
  filter(host_species == "Rattus rattus")
asv_names <- tibble(asv_ID = unique(data_asv_filtered$asv_ID))
asv_names %<>% mutate(ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", asv_ID))) %>% arrange(ID)
seq_fa2 <- seq_fa[names(seq_fa) %in% asv_names$asv_ID]
names(seq_fa2) <- asv_names$asv_ID
#write.FASTA(seq_fa2, file = "data/data_raw/data_microbiome/ASV_merged_full_0.01.fa") # writing the filtered fasta file

# how many ASVs in the tree?
nrow(asv_names)

# align the sequences
seq_aligned <- readDNAStringSet("data/data_raw/data_microbiome/ASV_merged_full_0.01.fa")
aligned <- DECIPHER::AlignSeqs(seq_aligned)
seq_aligned2 <- as.DNAbin(aligned)
seq_aligned <- msa::msaClustalOmega("data/data_raw/data_microbiome/ASV_merged_full_0.01.fa", type="dna",
                                    auto=F,
                                    cluster=100,
                                    dealign=F,
                                    order="input",
                                    useKimura=T)
seq_aligned2 <- as.DNAbin(seq_aligned)
rownames(seq_aligned2) <- asv_names$asv_ID

dnaphydatAll <- phyDat(seq_aligned2, type="DNA", levels=NULL)

mt2 <- modelTest(dnaphydatAll,
                 model="all", multicore = T, mc.cores=2)

# optimize the tree
# the function takes the best model (according to BIC) and optimizes the parameters
best_tree <- pml_bb(mt2) 

# save the optimized tree. this is a pml object. to get the tree call: fit$tree
#saveRDS(best_tree, file = "results/phylo_tree_0.01_2.rds")
#best_tree <- readRDS(file = "output/phylogenetic_tree/phylo_tree_0.01.rds")

# plot the tree
phylo_tree <- best_tree$tree
