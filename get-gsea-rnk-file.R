options(warn=1)

g <- read.delim("/mnt/projects/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.tsv")
#plot(density(log(apply(g[,10:21],1,var))))
#gf <- g[log(apply(g[,10:21],1,var))>3 & !is.na(g$hgnc_symbol) & g$hgnc_symbol != "",]
gf <- g[!is.na(g$pvalue) & g$pvalue<=0.01 & !is.na(g$hgnc_symbol) & g$hgnc_symbol != "" & !is.na(g$log2FoldChange),]
gf <- gf[order(gf$log2FoldChange,decreasing=T),c("hgnc_symbol", "log2FoldChange")]
write.table(gf, "/mnt/projects/chrisi/results/gsea/diff_exp_genes.rnk", col.names=F, row.names=F, sep="\t", quote=F)