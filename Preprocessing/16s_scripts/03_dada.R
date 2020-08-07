library(dada2)
setwd("trimmed/")

samples <- sort(scan("../inputs/samples", what="character"))

forward_reads <- paste0(samples, "_R1_trimmed.fq.gz")
reverse_reads <- paste0(samples, "_R2_trimmed.fq.gz")

filtered_forward_reads <- paste0("../filtered/", samples, "_R1.filtered.fq.gz")
filtered_reverse_reads <- paste0("../filtered/", samples, "_R2.filtered.fq.gz")

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, 
                              reverse_reads, filtered_reverse_reads, 
                              maxEE=c(2,2), rm.phix=TRUE, 
                              multithread=TRUE, minLen=150, 
                              truncLen=c(220, 180))


err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

dada_forward <- dada(derep_forward, err=err_forward_reads, 
                     multithread = TRUE, pool = TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, 
                     multithread=TRUE, pool = TRUE)


merged_amplicons <- mergePairs(dada_forward, derep_forward, 
                               dada_reverse, derep_reverse, 
                               trimOverhang=TRUE, minOverlap=50)

seqtab <- makeSequenceTable(merged_amplicons)

seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) 

saveRDS(seqtab_nochim, "../outputs/seqtab_nochim.rds")
