
# demultiplex with sabre --------------------------------------------------
# # to install sabre:
# # git clone https://github.com/najoshi/sabre.git
# # cd sabre
# # make

barcodes <- read.csv("inputs/16s_primer_barcodes.csv", stringsAsFactors = F)
barcodes <- barcodes[1:(nrow(barcodes) -1) , c("sequence.name", "blendID")]
colnames(barcodes) <- c("barcode", "r1")
barcodes$barcode <- gsub("GTGCCAGCMGCCGCGGTAA", "", barcodes$barcode)
barcodes$r1 <- paste0(barcodes$r1, "_R1.fq")
barcodes$r1 <- gsub("mock community #", "mockcomm", barcodes$r1)
barcodes$r1 <- gsub("negative control #", "negctrl", barcodes$r1)
barcodes$r1 <- gsub("negative control#", "negctrl", barcodes$r1)
barcodes$r2 <- gsub("_R1", "_R2", barcodes$r1)
write.table(barcodes, "inputs/16s_primer_barcodes_sabre.tsv", row.names = F, quote = F, sep = "\t", col.names = F)


# 
# # ~/bin/sabre/sabre pe -f 16S-16-17-Must_S1_L001_R1_001.fastq.gz -r 16S-16-17-Must_S1_L001_R2_001.fastq.gz -b ../16s_primer_barcodes_sabre.tsv -u no_bc_match_R1.fq -w no_bc_match_R2.fq
