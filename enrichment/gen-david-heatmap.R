options(error=recover)
source("~/chrisi/scripts/enrichment/pGSEAfunctions.R")
source("~/chrisi/scripts/enrichment/clusterFctCatDAVID.v3.R")

carps <- clusterFctCatDAVID.v3(
		Files=c("/home/STANNANET/christian.frech/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt",
				"/home/STANNANET/christian.frech/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt"),
		whichcols=c(12),
		whichPcols=c(12),
		PvalornumRow="numRow",
		Pval=0.05,
		howmany=100,
		largeCat=100000,
		Plot="NO",
		carps=NULL,
		pdftif="pdf",
		cut.col.scheme=0.4)

clusterFctCatDAVID.v3(
		Files=c("/home/STANNANET/christian.frech/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt",
				"/home/STANNANET/christian.frech/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt"),
		whichcols=c(12),
		whichPcols=c(12),
		PvalornumRow="numRow",
		Pval=0.05,
		howmany=100,
		largeCat=100000,
		Plot="YES",
		carps=carps,
		pdftif="pdf",
		cut.col.scheme=0.4)
