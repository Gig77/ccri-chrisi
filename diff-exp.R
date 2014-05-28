options(warn=1)

samples <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples[nrow(samples)+1,] <- c("18378", "18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count", "empty")
samples[nrow(samples)+1,] <- c("18379", "18379_AAGACA_C4993ACXX_5_20140509B_20140509.count", "empty")
samples[nrow(samples)+1,] <- c("18380", "18380_TAATCG_C4993ACXX_5_20140509B_20140509.count", "etv6")
samples[nrow(samples)+1,] <- c("18381", "18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count", "etv6")

library("biomaRt")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75

#--------
# DESeq
#--------
library("DESeq")

cds <- newCountDataSetFromHTSeqCount(samples, directory="htseq")
cds <- estimateSizeFactors(cds)
sizeFactors(cds)

cds <- estimateDispersions(cds)

res <- nbinomTest(cds, "empty", "etv6")

#write.table(res, file="strobl-dox-empty-vs-etv6.deseq.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#res <- read.delim(file="strobl-dox-empty-vs-etv6.deseq.tsv")
genes <- getGene(res$id, "ensembl_gene_id", mart)

res.annotated <- merge(res, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="strobl-dox-empty-vs-etv6.deseq.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#--------
# DESeq2
#--------
library("DESeq2")

cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory="htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)


genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="strobl-dox-empty-vs-etv6.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)
