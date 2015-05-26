
fuka <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
chrisi.strobl.empty_vs_etv6 <- read.delim(file="/mnt/projects/chrisi/results/strobl-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv")

m.strobl.empty_vs_etv6 <- merge(chrisi.strobl.empty_vs_etv6, fuka, by.x="hgnc_symbol", by.y="syms", all.x=T, all.y=T)
write.table(m.strobl.empty_vs_etv6, file="strobl-dox-empty-vs-etv6.deseq2.fuka.tsv", col.names=T, row.names=F, sep="\t", quote=F)

do_plot <- function(title) {
	fit <- lm(log2FoldChange~logFC, data=m.sig)
	p <- anova(fit)$'Pr(>F)'[1]
	R <- summary(fit)$r.squared
	
	plot(log2FoldChange~logFC, data=m.sig, xlab="log2FC Gerhard", ylab="log2FC Chrisi", ylim=c(-3, 3.5), xlim=c(-1.5, 3), main=sprintf("%s %.2g (R=%.2f, p=%.2g)", title, cutoff, R, p), cex=0.2)
	abline(fit, col="red")
	abline(v=0, lty=3)
	abline(h=0, lty=3)
	text(m.sig$logFC, m.sig$log2FoldChange-0.1, m.sig$hgnc_symbol, cex=0.5, col=(as.numeric(!is.na(m.sig$mausSym)) + as.numeric(!is.na(m.sig$Tijssen_Runx1_alone.1)) + as.numeric(!is.na(m.sig$Wi_hgnc.1)))+1)
}

pdf("/mnt/projects/chrisi/results/fuka-chrisi-scatter.significant-both.pdf")
for (cutoff in c(0.05, 0.01, 0.001, 0.0001)) {
	m.sig <- m.strobl.empty_vs_etv6[!is.na(m.strobl.empty_vs_etv6$Pval) & !is.na(m.strobl.empty_vs_etv6$pval) & m.strobl.empty_vs_etv6$Pval<=cutoff & m.strobl.empty_vs_etv6$pval<=cutoff,]
	do_plot("Both datasets p-value <=")
}
dev.off()

pdf("/mnt/projects/chrisi/results/fuka-chrisi-scatter.significant-chrisi.pdf")
for (cutoff in c(0.05, 0.01, 0.001, 0.0001, 0.00001, 1e-6, 1e-10)) {
	m.sig <- m.strobl.empty_vs_etv6[!is.na(m.strobl.empty_vs_etv6$pvalue) & !is.na(m.strobl.empty_vs_etv6$pval) & m.strobl.empty_vs_etv6$pvalue<=cutoff,]
	do_plot("Chrisi dataset p-value <=")
}
dev.off()