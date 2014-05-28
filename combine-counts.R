options(warn=1)

#args <- commandArgs(trailingOnly = TRUE)
args <- args <- c("~/chrisi/results/htseq/test.count", "~/chrisi/results/htseq/18376_ACATTA_C4993ACXX_5_20140509B_20140509.count")

if (is.na(args[1])) stop("ERROR: no count file specified")

c <- read.delim(args[1], header=F, colClasses=c("character", "integer"))
c <- c[!c$V1 %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"),]
colnames(c) <- c("gene", sub(".count$", "", basename(args[1])))

if (length(args) > 1) {
	for(i in 2:length(args)) {
		c1 <- read.delim(args[i], header=F, colClasses=c("character", "integer"))
		c1 <- c1[!c1$V1 %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"),]
		colnames(c1) <- c("gene", sub(".count$", "", basename(args[i])))
		c <- merge(c, c1, all.x=T, all.y=T)
	}
}

write.table(c, file=stdout(), col.names=T, row.names=F, sep="\t", quote=F)