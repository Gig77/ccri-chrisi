options(error=recover)
source("/mnt/projects/chrisi/scripts/enrichment/pGSEAfunctions.R")
source("/mnt/projects/chrisi/scripts/enrichment/clusterFctCatDAVID.v3.R")

carps <- clusterFctCatDAVID.v3(
		Files=c("/mnt/projects/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt",
				"/mnt/projects/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt"),
		whichcols=c(12),
		whichPcols=c(12),
		PvalornumRow="numRow",
		Pval=0.05,
		howmany=100,
		largeCat=100000,
		Plot="NO",
		carps=NULL,
		pdftif="pdf",
		cut.col.scheme=0.4,
		col.overRide="highsign")

clusterFctCatDAVID.v3(
		Files=c("/mnt/projects/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt",
				"/mnt/projects/chrisi/results/enrichment/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt"),
		whichcols=c(12),
		whichPcols=c(12),
		PvalornumRow="numRow",
		Pval=0.05,
		howmany=100,
		largeCat=100000,
		Plot="YES",
		carps=carps,
		pdftif="pdf",
		cut.col.scheme=0.4,
		col.overRide="highsign")
