
library(dada2)

seqtab_nochim <- readRDS("outputs/seqtab_nochim.rds")

taxa <- assignTaxonomy(seqtab_nochim, "inputs/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T)

# fix headers
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")

for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "outputs/ASVs.fa")

# count table:
asv_tab <- t(seqtab_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "outputs/ASVs_counts.txt", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "outputs/ASVs_taxonomy.txt", sep="\t", quote=F, col.names=NA)
