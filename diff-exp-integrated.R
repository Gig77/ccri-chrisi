options(warn=1)
library(optparse)
library("DESeq2")
library("biomaRt")

set.seed(343)

sample.table <- data.frame(name=character(), file=character(), expression=character(), construct=character(), vector=character(), dox=character(), stringsAsFactors=F)
sample.table[nrow(sample.table)+1,] <- c("18370", "18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count", "OFF", "ER", "zuber", "-")
sample.table[nrow(sample.table)+1,] <- c("18371", "18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count", "ON", "ER", "zuber", "+")
sample.table[nrow(sample.table)+1,] <- c("18372", "18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count", "OFF", "ER", "zuber", "-")
sample.table[nrow(sample.table)+1,] <- c("18373", "18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count", "ON", "ER", "zuber", "+")
sample.table[nrow(sample.table)+1,] <- c("18374", "18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.count", "OFF", "empty", "zuber", "-")
sample.table[nrow(sample.table)+1,] <- c("18375", "18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count", "OFF", "empty", "zuber", "+")
sample.table[nrow(sample.table)+1,] <- c("18376", "18376_ACATTA_C4993ACXX_5_20140509B_20140509.count", "OFF", "empty", "zuber", "-")
sample.table[nrow(sample.table)+1,] <- c("18377", "18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count", "OFF", "empty", "zuber", "+")
sample.table[nrow(sample.table)+1,] <- c("18378", "18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count", "OFF", "empty", "strobl", "+")
sample.table[nrow(sample.table)+1,] <- c("18379", "18379_AAGACA_C4993ACXX_5_20140509B_20140509.count", "OFF", "empty", "strobl", "+")
sample.table[nrow(sample.table)+1,] <- c("18380", "18380_TAATCG_C4993ACXX_5_20140509B_20140509.count", "ON", "ER", "strobl", "+")
sample.table[nrow(sample.table)+1,] <- c("18381", "18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count", "ON", "ER", "strobl", "+")

sample.table$expression <- factor(sample.table$expression, levels=c("OFF", "ON"))
sample.table$construct <- factor(sample.table$construct, levels=c("empty", "ER"))
sample.table$vector <- factor(sample.table$vector)
sample.table$dox <- factor(sample.table$dox, levels=c("-", "+"))

#----------------------------------------
# multi-factor
#----------------------------------------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table, directory="htseq", design= ~ construct + vector + dox + expression)
cds <- estimateSizeFactors(cds)
#sizeFactors(cds)
counts.norm <- as.data.frame(counts(cds, normalized=T))
dds <- DESeq(cds)
res <- results(dds)
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

# annotate genes with Ensembl biomart
#---
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "band", "strand", "start_position", "end_position"), mart=mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T) # add gene annotation
res.annotated <- res.annotated[,c(1,8,9,2,3,4,5,6,7)] # reorder columns
res.annotated <- merge(res.annotated, counts.norm, by.x="id", by.y="row.names", all.x=T)  # add normalized read counts to output
res.annotated <- res.annotated[order(res.annotated$padj),]

# write output
#---
write.table(res.annotated, file="integrated-expressed-vs-notexpressed.multi-factor.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#----------------------------------------
# only expression
#----------------------------------------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table, directory="htseq", design= ~ expression)
cds <- estimateSizeFactors(cds)
#sizeFactors(cds)
counts.norm <- as.data.frame(counts(cds, normalized=T))
dds <- DESeq(cds)
res <- results(dds, cooksCutoff=FALSE) # turn off filtering
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

# annotate genes with Ensembl biomart
#---
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "band", "strand", "start_position", "end_position"), mart=mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T) # add gene annotation
res.annotated <- res.annotated[,c(1,8,9,2,3,4,5,6,7)] # reorder columns
res.annotated <- merge(res.annotated, counts.norm, by.x="id", by.y="row.names", all.x=T)  # add normalized read counts to output
res.annotated <- res.annotated[order(res.annotated$padj),]

# write output
#---
write.table(res.annotated, file="integrated-expressed-vs-notexpressed.single-factor.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)
