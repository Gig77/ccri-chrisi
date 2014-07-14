options(warn=1)

samples.strobl <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.strobl[nrow(samples.strobl)+1,] <- c("18378", "18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count", "strobl-empty")
samples.strobl[nrow(samples.strobl)+1,] <- c("18379", "18379_AAGACA_C4993ACXX_5_20140509B_20140509.count", "strobl-empty")
samples.strobl[nrow(samples.strobl)+1,] <- c("18380", "18380_TAATCG_C4993ACXX_5_20140509B_20140509.count", "strobl-etv6")
samples.strobl[nrow(samples.strobl)+1,] <- c("18381", "18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count", "strobl-etv6")

samples.zuber <- data.frame(name=character(), file=character(), vector=character(), treatment=character(), stringsAsFactors=F)
samples.zuber[nrow(samples.zuber)+1,] <- c("18370", "18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count", "etv6", "nodox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18372", "18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count", "etv6", "nodox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18371", "18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18373", "18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18375", "18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count", "empty", "dox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18377", "18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count", "empty", "dox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18374", "18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.count", "empty", "nodox")
samples.zuber[nrow(samples.zuber)+1,] <- c("18376", "18376_ACATTA_C4993ACXX_5_20140509B_20140509.count", "empty", "nodox")

samples.zuber.etv6.doxVsNodox <- data.frame(name=character(), file=character(), vector=character(), treatment=character(), stringsAsFactors=F)
samples.zuber.etv6.doxVsNodox[nrow(samples.zuber.etv6.doxVsNodox)+1,] <- c("18370", "18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count", "etv6", "nodox")
samples.zuber.etv6.doxVsNodox[nrow(samples.zuber.etv6.doxVsNodox)+1,] <- c("18372", "18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count", "etv6", "nodox")
samples.zuber.etv6.doxVsNodox[nrow(samples.zuber.etv6.doxVsNodox)+1,] <- c("18371", "18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")
samples.zuber.etv6.doxVsNodox[nrow(samples.zuber.etv6.doxVsNodox)+1,] <- c("18373", "18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")

samples.zuber.dox.etv6VsEmpty <- data.frame(name=character(), file=character(), vector=character(), treatment=character(), stringsAsFactors=F)
samples.zuber.dox.etv6VsEmpty[nrow(samples.zuber.dox.etv6VsEmpty)+1,] <- c("18371", "18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")
samples.zuber.dox.etv6VsEmpty[nrow(samples.zuber.dox.etv6VsEmpty)+1,] <- c("18373", "18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count", "etv6", "dox")
samples.zuber.dox.etv6VsEmpty[nrow(samples.zuber.dox.etv6VsEmpty)+1,] <- c("18375", "18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count", "empty", "dox")
samples.zuber.dox.etv6VsEmpty[nrow(samples.zuber.dox.etv6VsEmpty)+1,] <- c("18377", "18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count", "empty", "dox")

library("biomaRt")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75

#--------
# DESeq2
#--------
library("DESeq2")

# STROBL ------------------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.strobl, directory="htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/chrisi/results/strobl-dox-empty-vs-etv6.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# ZUBER ETV6 DOX vs. NODOX ------------------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.zuber.etv6.doxVsNodox, directory="~/chrisi/results/htseq", design=~treatment)
dds <- DESeq(cds)
colData(dds)
dds$treatment <- relevel(dds$treatment, "nodox")
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/chrisi/results/zuber-etv6-dox-vs-nodox.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# ZUBER DOX ETV6 vs. EMPTY vector ------------------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.zuber.dox.etv6VsEmpty, directory="~/chrisi/results/htseq", design=~vector)
dds <- DESeq(cds)
colData(dds)
dds$vector <- relevel(dds$vector, "empty")
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/chrisi/results/zuber-dox-empty-vs-etv6.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#--------
# DESeq
#--------
library("DESeq")

cds <- newCountDataSetFromHTSeqCount(samples.strobl, directory="htseq")
cds <- estimateSizeFactors(cds)
sizeFactors(cds)

cds <- estimateDispersions(cds)

res <- nbinomTest(cds, "strobl-empty", "strobl-etv6")

#write.table(res, file="strobl-dox-empty-vs-etv6.deseq.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#res <- read.delim(file="strobl-dox-empty-vs-etv6.deseq.tsv")
genes <- getGene(res$id, "ensembl_gene_id", mart)

res.annotated <- merge(res, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="strobl-dox-empty-vs-etv6.deseq.tsv", col.names=T, row.names=F, sep="\t", quote=F)
