
#########################################################
cd ~/DATA/amnesia/Eleni/MouseRNASeq

mkdir FASTQC
fastqc   -t 6  -o FASTQC  *.fastq.gz  # outputdir hier schon angegeben
#########################################################



# cutadapt and fastx examples ###########################
cutadapt $(<cutadapt.conf) wt_David.fastq > Logwt.txt   2>&1 >wt_David.cuta.fastq   &

### fastx_trimmer
# fastx_trimmer example: trim first 5 bases
fastx_trimmer -f 5 -Q 33 -i dox_48h_David.cuta.fastq -o dox_48h_David.cuta1.fastq

# fastq_quality_filter example:
fastq_quality_filter -Q33 -i wt_David.cutaT.fastq  -o wt_David.cutaQT.fastq  -p 80  -q 25 -v  &


## example trim to 26 bases (like miRanalyzer)
fastx_trimmer -Q 33 -l 26  -i test.fastq  -o test26.fastq
##########################################################



### script including fastx-tools, gsnap, sort bam file ###
## align.sh ##############################################
#!/bin/sh
for i in *.bam  
   do
     echo $i
     Name=`echo "${i}" | sed 's/.unmapped.bam//g'`
     echo $Name
     fastq=${Name}.fastq
     fastqT=${Name}_qtrim.fastq

     #done

     java -jar ~/DATA/NGS/picard-tools-1.95/SamToFastq.jar  INPUT=$i  FASTQ=$fastq
     echo $fastq
     fastx_trimmer -f 8 -m 22 -Q 33 -i $fastq  | fastq_quality_filter -Q33  -p 80  -q 25 -o $fastqT  -v
     echo $fastqT

     rm $fastq
     #ohne mismatch - counted aber SNPs nicht als mismatch
     #out=${Name}.GSNAPOUT_nomism.sam
     #gsnap -A sam --batch=4  -n 5 -m 0 --quality-protocol=sanger --print-snps  -Q  -O  --nofails  -t 5  --input-buffer-size=3000  -d  hg19  -D /home/STANNANET/maximilian.kauer/DATA/NGS/GSNAP/hg19   -s hg19splicesitesfile   -v snpfile135  $i  >  $out  #2>&1 > gsnaplog.txt

     out=${Name}_GSNAPOUT.sam
     echo $out
     outb=${Name}_GSNAPOUT.bam
     echo $outb

     Sort=${Name}_GSNAPOUT_sort
     outbs=${Name}_GSNAPOUT_sort.bam
     #done

     #gsnap  -A sam --batch=4  --quality-protocol=sanger --print-snps  -Q --nofails  -t 5  --input-buffer-size=5000  -d  mm10  -D /home/STANNANET/maximilian.kauer/DATA/NGS/GSNAP/mm10   -s mm10splicesitesfile   -v snpfile137  $fastqT  >  $out
     gsnap  -A sam  -m 1 -n 1 -Q --nofails --batch=4  --quality-protocol=sanger --print-snps  -t 5  --input-buffer-size=5000  -d  mm10  -D /home/STANNANET/maximilian.kauer/DATA/NGS/GSNAP/mm10   -s mm10splicesitesfile   -v snpfile137  $fastqT  >  $out
     rm $fastqT

     samtools view -Sb  $out  >  $outb
     rm  $out     

     samtools sort $outb  $Sort
     samtools index $outbs
     rm $outb

     ## make coverage plots ###
     outt=${Name}.tdf
     outw=${Name}.wig
   
     echo $outt
     igvtools count -z 7 -w 15 -e 0 --strands read  $i  $outt,outw  mm10          
done
########################################################




### differential expression ############################
R
pathA <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/Eleni/RNASeq/Analysis"
bamdir <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/Eleni/RNASeq/fastq"

pathmirAnnot <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/genome/AnnotationData/miR/mirBase"
pathFct <-  "/home/STANNANET/maximilian.kauer/DATA/amnesia/R.Functions.etc"

load(paste(pathFct, "SLmisc", sep="/"))
load(paste(pathFct, "max.func.env", sep="/"))

source( paste(pathFct, "load.packages.R", sep="/") )

Name <- "EleniRNASeq"
#library(easyRNASeq)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(Rsamtools)
library(mirbase.db)
library(rtracklayer)

ensembl = useMart("ensembl")
ensembl = useDataset( "hsapiens_gene_ensembl", mart=ensembl )
#write.table( listAttributes(ensembl), file="ensemblAttributes.xls",  col.names=T, row.names=F, sep="\t", quote=F )
att <- listAttributes( mart=ensembl ); head( att, 20 )


## Fct um 1 Line Bm-file zu machen ###
shrinkFctBm <-  function(x) {
                     if(dim(x)[1] > 1) {
                          y <- x[1,]
                          for (i in length(dim(x)[1]) ) {
                              if( length( unique( x[,i][ -which(x[,i] == "") ] ) >1) ) {
                                    y[,i] <- paste( unique( x[,i] ), collapse="|" )
                                    }
                              if( length( unique( x[,i][ -which(x[,i] == "") ] ) == 1) ) {
                                    y[,i] <- unique( x[,i][ -which(x[,i] == "") ] )
                                    }
                          }
                      }
                      if(dim(x)[1] == 1) {  y <- x  }
                      y
            }
########################################################


setwd(pathA)


## get annotation ######################################
# form biomart
hse <- makeTranscriptDbFromBiomart( biomart="ensembl",
            dataset="hsapiens_gene_ensembl", miRBaseBuild="GRCh37.p5" )
#!! geht auch:  makeTranscriptDbFromGFF
#cdsByGene <- cdsBy( hse, by="gene" ) ## git auch cdsBy, intronsBy, fiveUTRsByTranscript, threeUTRsByTranscript
exonsByGene <- exonsBy(hse, by="gene") ## git auch cdsBy, intronsBy, fiveUTRsByTranscript, threeUTRsByTranscript
#seqnames(cdsByGene)

newLev <- paste( "chr", seqlevels( exonsByGene )[1:25], sep="" )
names(newLev) <-  seqlevels( exonsByGene )[1:25]
exonsByGene1 <- renameSeqlevels( exonsByGene , newLev )
seqlevels( exonsByGene1 )[1:24]
exonsByGene2 <- keepSeqlevels( exonsByGene1 , seqlevels( exonsByGene1 )[1:24] )
seqlevels( exonsByGene2 )

save( exonsByGene2, file="exByGene.RObj" )
save( exonsByGene, file="exByGene_orig.RObj" )

saveDb( hse, file="txDbEnsembl.sqlite" )

# test
#setwd("test")
#rGA1 <- readGAlignmentsList( file=fls[1] )
########################################################



### summarize Overlaps #################################
setwd( pathA )
load( paste(bamdir, "exByGene.RObj",sep="/") ); exonsByGene <- exonsByGene2; rm( exonsByGene2 )
seqlevels( exonsByGene )
class( exonsByGene )

##
#setwd(bamdir)

pat <- ".bam$"
fls <- list.files( path=bamdir, pattern=pat, full=TRUE ) ; fls
bamlst <- BamFileList( fls, yieldSize=100000 ); bamlst

se <- summarizeOverlaps( exonsByGene, bamlst, mode="Union",
                       singleEnd=T, ignore.strand=T )

colData(se); colnames(se)
colData(se)$Name <- paste("Eleni", gsub( "_GSNAPOUT_sort.bam", "", sapply( strsplit( names(bamlst), "/" ), function(x) i=tail(x,1) ) ) , sep="")
colData(se)$treat <- "nodox"; colData(se)$treat[ grep("dox", colData(se)$Name) ] <- "dox"
colData(se)$treat

setwd( pathA )
save(se, file=paste( "summarizedOverlaps" , Name, sep="_") )

load( "summarizedOverlaps_EleniRNASeq" )
counts <-  assays(se)$counts; colnames(counts) <- colData(se)$Name
head(counts); apply( counts, 2, sum )
########################################################




## DEseq2 ##############################################
setwd(pathA)
library( DESeq2 )

dds <- DESeqDataSet(se = se, design = ~as.factor(treat) )

#dds <- DESeq(dds); CooksDist <- "cooksDistYES"  # problem wegen pvalues NA:
#The Cook's cutoff used here is the .75 quantile of the F(p, m-p) distribution. So for 2 parameters (intercept and condition) and 10 samples, we use a cutoff of:
#> qf(.75, 2, 10 - 2)
#[1] 1.656854
#The sample with count 60 is getting 2.211495, a high estimate of Cook's distance (meaning the log fold change would change a lot if this sample were removed),
#although from looking at the counts I agree that this is not desirable behavior to give this gene a p-value of NA.
#I will investigate how we might better set this cutoff.
#You can set a higher cutoff manually or turn off the filtering by Cook's with:

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# turn off filtering:
dds <- nbinomWaldTest(dds, cooksCutoff=FALSE) ; CooksDist <- "cooksDistNO"
res <- results(dds)

# or use a higher cutoff:
#p <- 2; m <- 4
#dds <- nbinomWaldTest( dds, cooksCutoff=qf(.95,p,m-p) )

# check if p values are there
mcols(dds)[171,]
res[171,]


dds <- estimateSizeFactors( dds )
sizeFactors( dds )

pdf( file=paste("MAplot",  "pdf", sep="."), height=9, width=12 )
  plotMA(dds)
graphics.off()

counts <-  assays(se)$counts; colnames(counts) <- colnames(counts) <- colData(se)$Name;
counts_norm <- counts
for( i in 1:ncol(counts_norm) ) {
       counts_norm[, i] <- round( counts_norm[ ,i]/sizeFactors(dds)[i] )
}
colnames(counts_norm) <- paste("norm",colnames(counts_norm),sep="")

stopifnot( identical( rownames(res), rownames(counts) ) )
matAnn <- cbind( as.data.frame( res, stringsAsFactors=F), as.data.frame( counts, stringsAsFactors=F), counts_norm )
matAnn$ensg <- rownames(matAnn); matAnn <- matAnn[ , c( ncol(matAnn), 1:(ncol(matAnn)-1) )]; head(matAnn)

write.table(matAnn, file=paste("resMat", Name, CooksDist, "xls", sep="."), quote=F, row.names=F, sep="\t")
###################################################################



## RPKM ###########################################################
geneLengths <- sum( end( exonsByGene ) - start( exonsByGene ) )/1000; head(geneLengths)
sumCounts <- apply( counts, 2, sum )/1000000 ; sumCounts

rpm <- as.data.frame( counts , stringsAsFactors=F )
for( i in 1:ncol(counts)) {
            rpm[, i] <- as.numeric( counts[, i]/sumCounts[i] )
}
head(rpm)

tmp <- merge( rpm, as.data.frame( geneLengths, stringsAsFactors=F ), by="row.names" ); head(tmp)
rpkm <- as.data.frame( as.matrix( tmp[, 2:5] ) / as.numeric( tmp[,6] ) )
rownames(rpkm) <- tmp$Row.names; colnames(rpkm) <-  paste("rpkm", colnames(rpkm), sep=""); head( rpkm )
###################################################################



## write counts out ###############################################
#out <- counts; outName <- "Counts"
setwd(pathA)
out <- as.data.frame( res ) ; outName <- paste("DESeqRes", Name, sep="_")
out <- out[ -which(out$baseMean == 0 ) ,] ; dim(out)

bm1 <- getBM(attributes=c("hgnc_symbol",  "ensembl_gene_id"), filters = "ensembl_gene_id", values=as.character(rownames(out) ), mart = ensembl); head(bm1)
bm2 <- getBM(attributes=c("rfam","refseq_ncrna","ensembl_gene_id"), filters = "ensembl_gene_id", values=as.character(rownames(out) ), mart = ensembl); head(bm2)
bm3 <- merge(bm1, bm2, by="ensembl_gene_id" ); bm3 <- bm3[order(bm3$ensembl_gene_id), ]; head(bm3);

whd <- sort( union( which( duplicated(bm3$ensembl_gene_id) ), which( duplicated(bm3$ensembl_gene_id, fromLast=T) ) ) ); head(whd)
bm <- bm3[ -whd, ]
system.time( test <- by( bm3[whd,], as.character(bm3[whd,1]), function(x) { x <- shrinkFctBm(x) } ) )
system.time(  df1 <- do.call( rbind, test )  )
bm <- rbind(bm, df1); head(bm); dim(bm)

out1 <- merge( bm, out,  by.x="ensembl_gene_id" , by.y="row.names",  all.x=F, all.y=T ); head(out1)
if( outName == "Counts") { out1 <- out1[ order(out1$IN_1, decreasing=T),] }
if( unlist(strsplit(outName,"_"))[1] == "DESeqRes") {
           mout1 <- merge( as.data.frame( out1 ), as.data.frame( counts_norm ), by.x="ensembl_gene_id", by.y="row.names" )
           mout1 <- merge( mout1, as.data.frame( counts ), by.x="ensembl_gene_id", by.y="row.names" )
           mout1 <- merge( mout1, as.data.frame( rpkm ), by.x="ensembl_gene_id", by.y="row.names" )
           mout1 <- mout1[ order(mout1$padj), ]
}

head(mout1); dim(mout1)

## make means
mout1$meanNormCountNodox <- rowMeans( mout1[, c("normEleniE29_1", "normEleniE31_3")] )
mout1$meanNormCountDox <- rowMeans( mout1[, c("normEleniE30_2_dox", "normEleniE32_4_dox")] )
mout1$meanRPKMNodox <- rowMeans( mout1[, c("rpkmEleniE29_1", "rpkmEleniE31_3")] )
mout1$meanRPKMDox <- rowMeans( mout1[, c("rpkmEleniE30_2_dox", "rpkmEleniE32_4_dox")] )
head(mout1); dim(mout1)

write.table( as.data.frame( mout1 ), file=paste(outName, CooksDist, "xls", sep="."), col.names=T, row.names=F, sep="\t", quote=F )
save( mout1, file=paste(outName, CooksDist, "RObj", sep=".") )
rm( bm1, bm2, bm3, bm, test, df1, out )
###################################################################



