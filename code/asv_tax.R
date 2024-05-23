# script for assignung taxonomy to ASVs

#BiocManager::install("dada2")
library(tidyverse)
library(dada2)
library(DECIPHER)
library(ape)
library(adegenet)
library(phangorn)

# reading the fasta file
fasta_file <- read.FASTA(file="data/data_raw/data_microbiome/ASV_merged_full.fa")


new_taxa = assignTaxonomy("data/data_raw/data_microbiome/ASV_merged_full.fa",
                          "~/Research/gut_parasites/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                          verbose=T, minBoot=80, outputBootstraps = T)


# Save the taxonomic output

# check the structure
str(new_taxa)

# number of asvs
n_asvs = dim(new_taxa[[1]])
new_ids = t(as.data.frame(new_taxa[[1]][1,1:7])) # check that there are 7 taxonomic slots

# loop through (this is not the most efficient way to do things, but it only takes a minute)
for(i in 2:n_asvs){
  temp = t(as.data.frame(new_taxa[[1]][i,1:7]))
  new_ids= rbind(new_ids, temp)
}
dim(new_ids)
new_ids = as.data.frame(new_ids)
names(rownames(new_taxa$tax))
# double check names
new_ids$ASV = names(rownames(new_taxa$tax))
# write out the file
write.table(new_ids, "data/data_raw/data_microbiome/ASVs_taxonomy_new.tsv", sep = "\t", quote=F, col.names=NA)


### Repeat this procedure for the boot scores
### Boot values
n_asvs = dim(new_taxa[[2]])
new_ids = t(as.data.frame(new_taxa[[2]][1,1:7]))

for(i in 2:n_asvs){
  temp = t(as.data.frame(new_taxa[[2]][i,1:7]))
  new_ids= rbind(new_ids, temp)
}
dim(new_ids)
str(new_ids)
new_ids = as.data.frame(new_ids)
new_ids$ASV = names(rownames(new_taxa$tax))


write.table(new_ids, "data/data_raw/data_microbiome/ASVs_taxonomy_new_boot.tsv", sep = "\t", quote=F, col.names=NA)
