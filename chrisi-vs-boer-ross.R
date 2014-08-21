boer <- read.delim("~/chrisi/data/RossBoer/NordischALL.esetnsF.annot.txt")
boer.2009 <- read.delim("~/chrisi/data/RossBoer/matAnn.GSE13351_BOER.eset_zfilt_th3_nsF.tsv")
ross <- read.delim("~/chrisi/data/RossBoer/ROSS2.2003.esetnsF.annot.txt")
chrisi <- read.delim("~/chrisi/results/zuber-dox-empty-vs-etv6.deseq2.tsv")

m.ross <- merge(ross, chrisi, by.x="syms", by.y="hgnc_symbol")
m.boer <- merge(boer, chrisi, by.x="syms", by.y="hgnc_symbol")
m.boer.2009 <- merge(boer.2009, chrisi, by.x="syms", by.y="hgnc_symbol")

gene.highlights <- c("ETV6", "RUNX1", "EPOR")

pdf("~/chrisi/results/chrisi-boer-scatter.pdf")
fit <- lm(log2FoldChange~TAvs.mean.noTall, data=m.boer)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m.boer$TAvs.mean.noTall, m.boer$log2FoldChange, xlim=c(-4, 5), ylim=c(-4, 4), xlab="Boer (TAvs.mean.noTall)", ylab="Zuber (log2FoldChange)", main=sprintf("Zuber-dox-etv6-vs-empty vs. Boer (R=%.2f, p=%.2g)", R, p), cex=0.3)
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- m.boer$adjPval.TAvs.mean.noTall <= 0.1 & m.boer$padj <= 0.1
sig[is.na(sig)] <- FALSE
text(m.boer$TAvs.mean.noTall[!sig], m.boer$log2FoldChange[!sig]-0.1, m.boer$syms[!sig], cex=0.3, col="black")
text(m.boer$TAvs.mean.noTall[sig], m.boer$log2FoldChange[sig]-0.1, m.boer$syms[sig], cex=0.3, col="red")
text(m.boer$TAvs.mean.noTall[m.boer$syms %in% gene.highlights], m.boer$log2FoldChange[m.boer$syms %in% gene.highlights]-0.1, m.boer$syms[m.boer$syms %in% gene.highlights], cex=0.3, col="orange")
dev.off()

pdf("~/chrisi/results/chrisi-boer2009-scatter.pdf")
fit <- lm(log2FoldChange~TA_vs_rest, data=m.boer.2009)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m.boer.2009$TA_vs_rest, m.boer.2009$log2FoldChange, xlim=c(-4, 5), ylim=c(-4, 4), xlab="Boer 2009 (TA_vs_rest)", ylab="Zuber (log2FoldChange)", main=sprintf("Zuber-dox-etv6-vs-empty vs. Boer 2009 (R=%.2f, p=%.2g)", R, p), cex=0.3)
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- m.boer.2009$adjP.TA_vs_rest <= 0.1 & m.boer.2009$padj <= 0.1
sig[is.na(sig)] <- FALSE
text(m.boer.2009$TA_vs_rest[!sig], m.boer.2009$log2FoldChange[!sig]-0.1, m.boer.2009$syms[!sig], cex=0.3, col="black")
text(m.boer.2009$TA_vs_rest[sig], m.boer.2009$log2FoldChange[sig]-0.1, m.boer.2009$syms[sig], cex=0.3, col="red")
text(m.boer.2009$TA_vs_rest[m.boer.2009$syms %in% gene.highlights], m.boer.2009$log2FoldChange[m.boer.2009$syms %in% gene.highlights]-0.1, m.boer.2009$syms[m.boer.2009$syms %in% gene.highlights], cex=0.3, col="orange")
dev.off()

pdf("~/chrisi/results/chrisi-ross-scatter.pdf")
fit <- lm(log2FoldChange~TAvs.mean_noTALL, data=m.ross)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m.ross$TAvs.mean_noTALL, m.ross$log2FoldChange, xlim=c(-4, 5), ylim=c(-4, 4), xlab="Ross (TAvs.mean_noTALL)", ylab="Zuber (log2FoldChange)", main=sprintf("Zuber-dox-etv6-vs-empty vs. Ross (R=%.2f, p=%.2g)", R, p), cex=0.3)
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- m.ross$adjPval.TAvs.mean_noTALL <= 0.1 & m.ross$padj <= 0.1
sig[is.na(sig)] <- FALSE
text(m.ross$TAvs.mean_noTALL[!sig], m.ross$log2FoldChange[!sig]-0.1, m.ross$syms[!sig], cex=0.3, col="black")
text(m.ross$TAvs.mean_noTALL[sig], m.ross$log2FoldChange[sig]-0.1, m.ross$syms[sig], cex=0.3, col="red")
text(m.ross$TAvs.mean_noTALL[m.ross$syms %in% gene.highlights], m.ross$log2FoldChange[m.ross$syms %in% gene.highlights]-0.1, m.ross$syms[m.ross$syms %in% gene.highlights], cex=0.3, col="orange")
dev.off()
