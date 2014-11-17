chrisi <- read.delim("~/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.tsv", check.names=F)

boer <- read.delim("~/chrisi/data/RossBoer/NordischALL.esetnsF.annot.txt", check.names=F)
boer.2009 <- read.delim("~/chrisi/data/RossBoer/matAnn.GSE13351_BOER.eset_zfilt_th3_nsF.tsv", check.names=F)
ross <- read.delim("~/chrisi/data/RossBoer/ROSS2.2003.esetnsF.annot.txt", check.names=F)
fuka <- read.delim("~/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv", check.names=F)

boer.short <- boer[,c("syms", "TAvs.mean.noTall", "adjPval.TAvs.mean.noTall")]
names(boer.short) <- c("syms", "boer.TAvs.mean.noTall", "boer.adjPval.TAvs.mean.noTall")
boer.2009.short <- boer.2009[,c("syms", "TA_vs_rest", "adjP.TA_vs_rest")]
names(boer.2009.short) <- c("syms", "boer2009.TA_vs_rest", "boer2009.adjP.TA_vs_rest")
ross.short <- ross[,c("syms", "TAvs.mean_noTALL", "adjPval.TAvs.mean_noTALL")]
names(ross.short) <- c("syms", "ross.TAvs.mean_noTALL", "ross.adjPval.TAvs.mean_noTALL")
fuka.short <- fuka[,c("syms", "logFC", "Padj")]
names(fuka.short) <- c("syms", "fuka.logFC", "fuka.Padj")

m <- merge(chrisi, fuka.short, by.x="hgnc_symbol", by.y="syms", all.x=T, all.y=T)
m <- merge(m, boer.short, by.x="hgnc_symbol", by.y="syms", all.x=T, all.y=T)
m <- merge(m, boer.2009.short, by.x="hgnc_symbol", by.y="syms", all.x=T, all.y=T)
m <- merge(m, ross.short, by.x="hgnc_symbol", by.y="syms", all.x=T, all.y=T)

write.table(m, file="~/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.fuka-ross-boer.tsv", col.names=T, row.names=F, sep="\t", quote=F)
