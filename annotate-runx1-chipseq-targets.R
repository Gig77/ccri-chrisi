


#### vgl mit ChIpseq ############################################################
########## Tijssen-Gottwald ################################
## read komische Tabelle
setwd(paste(pathBal, "ChIPseqOverlap", sep="/"))
tj <-  read.csv("chipseq/Tijssen_all.genes.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T); tj[1:5,]
colnames(tj) <-  paste("Tijssen",colnames(tj),sep="_")
allRunx1 <- as.vector(unique(do.call(c, tj)))


########## Wilson-Gottgens #################################
wi <-  read.csv("chipseq/Wilson_Gottgens_ChIPseq.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T); wi[1:5,]
wiRunx1 <- unique(wi$Runx1); wiRunx1 <- wiRunx1[ -which(wiRunx1 == "") ]; wiRunx1

## get orthologs
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
humOrt <- getLDS(attributes = c("mgi_symbol"),
       filters = "mgi_symbol", values=wiRunx1, mart = mouse,
       attributesL = c("hgnc_symbol", "entrezgene"), martL = human )
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[-which(is.na(humOrt$entrezgene)),]
wihumOrt <- humOrt; colnames(wihumOrt) <- c("Wi_mgi_symbol", "Wi_hgnc", "Wi_entrezgene")
wiRunx1hu <- sort( unique(humOrt$hgnc[humOrt$hgnc != ""]) ); wiRunx1hu
#############################################################


#if(interactive()){
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"), filters = "hgnc_symbol", values = "TP53", mart = human, attributesL = c("chromosome_name","start_position"), martL = mouse)
#}


########## Niebuhr ##########################################
ni <-  read.csv("chipseq/Niebuhr_TableS3_Runx1 Peaks Called in ProB-Cells.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T); ni[1:5,]
par(mfrow=c(2,2)); hist(ni$dist_tss);  hist(ni$score); plot(ni$dist_tss,ni$score, pch=20, cex=0.8); plot(ni$dist_tss,ni$score, xlim=c(-50000, 50000), pch=20, cex=0.8)


cutoffDist <- c(-5000, 1000)
nif <- ni[ which(ni$dist_tss > cutoffDist[1] & ni$dist_tss < cutoffDist[2]), ]; hist(nif$dist); hist(nif$score, br=100); length( unique( nif$nearest.gene))
plot(nif$dist_tss,nif$score, pch=20, cex=0.8);

scoreCutoff <- 100
nifs <- nif[which(nif$score > scoreCutoff)  ,]; hist(nifs$dist); hist(nifs$score, br=100); length( unique( nifs$nearest.gene))

mausIDs <- unique(ni$nearest.gene); tName <- "allNiebuhr"
#mausIDs <- nif$nearest.gene; tName <- "noScorecutoff"
#mausIDs <- nifs$nearest.gene; tName <- scoreCutoff

## get orthologs
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
humOrt <- getLDS(attributes = c("mgi_symbol"),
       filters = "mgi_symbol", values=mausIDs, mart = mouse,
       attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene"); head(humOrt)
if( length(which(is.na(humOrt$entrezgene)))>0 ) { humOrt <- humOrt[-which(is.na(humOrt$entrezgene)),] }
rm(mausIDs)

if(tName == "allNiebuhr") { allnieRunx1hu <- sort(unique(humOrt$hgnc[humOrt$hgnc != ""])); allnieRunx1hu }
if(tName == "noScorecutoff") { nieRunx1hu <- sort(unique(humOrt$hgnc[humOrt$hgnc != ""])); nieRunx1hu }
if(tName == scoreCutoff ) { nieScore100Runx1hu <- sort(unique(humOrt$hgnc[humOrt$hgnc != ""])); nieScore100Runx1hu }

uniqOrt <- unique(humOrt[,1:2])
By <- by( uniqOrt, uniqOrt$mgi_symbol, function(x) { y <-  paste(sort(x$hgnc),collapse="|"); y }, simplify=F )
dfOrt <- as.data.frame(do.call(rbind, By))
dfOrt$V1 <- gsub("^\\|", "", dfOrt$V1)

By <- by( ni, ni$nearest.gene, function(x) { y <-  paste(sort(x$dist_tss),collapse="|"); y }, simplify=F )
dfni <- as.data.frame(do.call(rbind, By)); head(dfni)

ni1Line <- merge(dfni, dfOrt, by="row.names", all=F); head(ni1Line)
colnames(ni1Line) <- c("mausSym", "Niedist","menschSym")

#write.table(ni1Line , file=paste(paste("MausHuman", "Niebuhr1Line_all", "xls", sep="."), sep="/"), quote=F, row.names=F, col.names=F, na="", sep="\t")
#############################################################



## ChIpseq an matAnn dranh�ngen #############################
chrisi <- read.delim(file="~/chrisi/results/strobl-dox-empty-vs-etv6.deseq2.tsv", na.strings=c("NA", ""))
chrisi$hgnc_symbol[chrisi$hgnc_symbol==""] <- NA
merged <- merge(chrisi, ni1Line, by.x="hgnc_symbol", by.y="menschSym", all.x=T, all.y=F); head(merged); dim(merged)
merged <- merge(merged, tj[,c(3,3)], by.x="hgnc_symbol", by.y="Tijssen_Runx1_alone", all.x=T, all.y=F); head(merged); dim(merged)
wihumOrt.dedup <- wihumOrt[,c(2,2)]
wihumOrt.dedup <- wihumOrt.dedup[!duplicated(wihumOrt.dedup),]
merged <- merge(merged, wihumOrt.dedup, by.x="hgnc_symbol", by.y="Wi_hgnc", all.x=T, all.y=F); head(merged); dim(merged)

merged <- merged[order(merged$padj), ]
write.table(merged, file="~/chrisi/results/strobl-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv", col.names=T, row.names=F, sep="\t", quote=F)
#############################################################

