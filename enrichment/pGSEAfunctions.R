########################### FUNCTIONS ################################################
######################################################################################
######################################################################################

######
#library(PGSEA)
#pathGSEA <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/genome/GSEA.GMTs"
#c1 <- readGmt(paste(pathGSEA,"c1.all.v2.5.symbols.gmt",sep="/") )
#c2 <- readGmt(paste(pathGSEA,"c2.all.v2.5.symbols.gmt",sep="/") )
#c2p <- readGmt(paste(pathGSEA,"c2.cp.v2.5.symbols.gmt",sep="/") )
#c3 <- readGmt(paste(pathGSEA,"c3.all.v2.5.symbols.gmt",sep="/") )
#c3tf <- readGmt(paste(pathGSEA,"c3.tft.v2.5.symbols.gmt",sep="/") )
#c4 <- readGmt(paste(pathGSEA,"c4.all.v2.5.symbols.gmt",sep="/") )
#c5 <- readGmt(paste(pathGSEA,"c5.all.v2.5.symbols.gmt",sep="/") )
#gmtlist <- list(c1,c2,c2p,c3,c3tf, c4, c5); names(gmtlist) <- c("positional gene sets","curatedgenesets","pathways","motifgenesets","TFmotifgenesets", "computationalgenesets", "GO")
#namesshort <- c("c1","c2","c2p","c3","c3tf","c4","c5")
######
#######################################################################################



####################################################
makePGSEAmat <- function(mat, colNames) {

    grb <- mat[,c(colNames)];
    colnames(grb)[1] <- "syms"; head(grb)

    whdup <- unique(c(which(duplicated(grb$syms)), which(is.na(grb$syms))));
    if (length(whdup)>0) { grbsub <- grb[-whdup,] }
    if (length(whdup)==0) { grbsub <- grb }
    rownames(grbsub) <- grbsub$syms; head(grbsub)

    pgseamat <- grbsub; pgseamat <- as.data.frame(pgseamat[,-1],stringsAsFactors=F);
    head(pgseamat); dim(pgseamat)
    pgseamat
}
####################################################





#### makeDF ########################################
makeDF <- function(pslist) {
    DF <- data.frame(GO="GO", Genes="Genes", stringsAsFactors=FALSE)
    for (i in 1:length(pslist)) {
         if (length(pslist[[i]]) > 1) {
         df <- data.frame(GO=names(pslist)[[i]], Genes=as.vector(pslist)[[i]])
         DF <- rbind(DF,df)
         }
    }; DF <- DF[-1,]; #DF[1:5,]
return (DF)
}
####################################################


####################################################
writeGSEA <- function(PG, pth, intersection) {
    for (i in 1:length(PG)) {
          pgdf <- PG[[i]];
          howmany <- dim(pgdf)[[1]]
          pgdf$Name <- rownames(pgdf); pgdf <- pgdf[,c(ncol(pgdf), 1:(ncol(pgdf)-1) )]
          #### an pg Genes dranhaengen zum ausdrucken f. TabS3
          pgdf$genes <- vector(length=dim(pgdf)[[1]])
          cl <- gmtlist[[i]];  names(cl) <- gsub(' na',"",names(cl));  names(cl) <- gsub('na',"",names(cl))
                               names(cl) <- gsub(" .+","",names(cl))
          for(j in 1:howmany) {
                rn <- rownames(pgdf)[j]
                isrn <- intersect(rownames(pgseamat), cl[[rn]]@ids);
                pgdf$genes[j] <- paste(isrn, collapse=",")
          }
       pgdf <- pgdf[order(pgdf[,2]),]
       colnames(pgdf) <- gsub("results\\.","",colnames(pgdf))
       write.table(pgdf, file=paste("GSEAtests", tabName, names(PG)[i],"xls",sep="."), sep="\t", quote=F, row.names=F)
    
        ## filter pgdf for all columns sign. ##
        if(intersection == "YES") {
              tmp <- pgdf[,c(-1,-ncol(pgdf))]
              whPsig <- which(as.vector(apply(tmp[,c((ncol(tmp)/2+1):ncol(tmp))],1, function(x) all(x < pth))))
              whCons <- which(as.vector(apply(tmp[,c(1:(ncol(tmp)/2))],1, function(x) sum(sign(range(x))) != 0))) 
              whBoth <- intersect(whPsig, whCons)
              if(length(whBoth)>0) { 
                  pgdfsig <- pgdf[whBoth,]
                  write.table(pgdfsig, file=paste("Bothsig.GSEAtests",pth, tabName, names(PG)[i],"xls",sep="."), sep="\t", quote=F, row.names=F)
                 }
        }
        #######################################    
    }
}
####################################################



####################################################
writeGSEAinList <- function(PG, pth, intersection) {
    outList <- list()
    for (i in 1:length(PG)) {
          pgdf <- PG[[i]];
          howmany <- dim(pgdf)[[1]]
          pgdf$Name <- rownames(pgdf); pgdf <- pgdf[,c(ncol(pgdf), 1:(ncol(pgdf)-1) )]
          #### an pg Genes dranhaengen zum ausdrucken f. TabS3
          pgdf$genes <- vector(length=dim(pgdf)[[1]])
          cl <- gmtlist[[i]];  names(cl) <- gsub(' na',"",names(cl));  names(cl) <- gsub('na',"",names(cl))
                               names(cl) <- gsub(" .+","",names(cl))
          for(j in 1:howmany) {
                rn <- rownames(pgdf)[j]
                isrn <- intersect(rownames(pgseamat), cl[[rn]]@ids);
                pgdf$genes[j] <- paste(isrn, collapse=",")
          }
       pgdf <- pgdf[order(pgdf[,2]),]
       colnames(pgdf) <- gsub("results\\.","",colnames(pgdf))
       outList[[i]] <- pgdf; names(outList)[i] <- names( PG )[i]        
       #write.table(pgdf, file=paste("GSEAtests", tabName, names(PG)[i],"xls",sep="."), sep="\t", quote=F, row.names=F)
    }
    outList
}
####################################################



#### PGSEA #########################################
pGSEAfun <-  function(pgseamat, Pval, weight, nRow, cRow,tName,absGSEA, groupSize) {

     #library(PGSEA)
     #pathGSEA <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/genome/GSEA.GMTs"
     #c1 <- readGmt(paste(pathGSEA,"c1.all.v2.5.symbols.gmt",sep="/") )
     #c2 <- readGmt(paste(pathGSEA,"c2.all.v2.5.symbols.gmt",sep="/") )
     #c3 <- readGmt(paste(pathGSEA,"c3.all.v2.5.symbols.gmt",sep="/") )
     #c4 <- readGmt(paste(pathGSEA,"c4.all.v2.5.symbols.gmt",sep="/") )
     #c5 <- readGmt(paste(pathGSEA,"c5.all.v2.5.symbols.gmt",sep="/") )
     #gmtlist <- list(c1,c2,c3, c4, c5); names(gmtlist) <- c("positional gene sets","curated gene sets ","motif gene sets", "computational gene sets", "GO")
     #namesshort <- c("c1","c2","c3","c4","c5")

    Pval <- Pval
    pgList <- list()
    for (i in 1:length(gmtlist)) {
        #i <- 3;  Pval=T; weight=F; nRow=50; cRow=0.4; tName=tabName; absGSEA="no"
        nSample <- ncol(pgseamat)
        nRow=nRow
        pg <- PGSEA(as.matrix(pgseamat), cl=gmtlist[[i]], range=groupSize, ref=NULL, center=FALSE, p.value=Pval, weighted=weight, enforceRange=TRUE)
        pg <- as.data.frame(pg)
        #if (Pval != T) {
                         pg[is.na(pg)] <- 0
                         pg <- pg[apply(pg,1,sum) != 0,] #}
        #pg <- pg[order(pg[,1]),]
        rownames(pg) <- gsub(' na',"",rownames(pg))  #geht nicht wegen overlap unten aber fuer smc-plots schon
        rownames(pg) <- gsub('na',"",rownames(pg))  #geht nicht wegen overlap unten aber fuer smc-plots schon
        rownames(pg) <- gsub(" .+","",rownames(pg))  
        
        head(pg)

        pdf(file=paste(names(gmtlist)[i],tName, "smcPlot",nRow,"mostsign","pdf", sep="."), height=12, width=8)
            #devSVG(file=paste(names(gmtlist)[i], Pval, "smcPlot.svg", sep="."), height=12, width=8)
            for(j in 1:nSample) {
            pg <- pg[ order(pg[,j]), ]

            if(dim(pg)[[1]] <= 2*nRow) { plotpg <- pg[,1:nSample] }
            if(dim(pg)[[1]] > 2*nRow) {
                  plotpg <- pg[ c(1:nRow,(nrow(pg)-nRow-1):nrow(pg)), 1:nSample ]
            }
            
            if(absGSEA == "no") { hcols <- SLmisc$heatmapCol(data=as.matrix(plotpg), col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(plotpg)))) ) }

            if(absGSEA == "yes") {
                       plotpg <- pg[c((nrow(pg)-(nRow*2)-1):nrow(pg)),1:nSample]
                       #hcols <- colorpanel(100, "yellowgreen","goldenrod", "tomato")
                       hcols <- SLmisc$heatmapCol(data=as.matrix(plotpg), col=colorpanel(100, "yellowgreen","orange", "tomato") , lim=min(abs(range(as.matrix(plotpg)))) )
            }
            #smcPlot(as.matrix(plotpg), scale=c(floor(min(plotpg)),ceiling(max(plotpg[,1:nSample]))), r.cex=0.5,c.cex=0.8, cex.main=0.6, main = paste(names(gmtlist)[i], colnames(pg)[j]), show.grid = T, margins = c(1,1, 8, 30), col=hcols)
            par(cex.main=0.5)
            heatmap.2(as.matrix(plotpg), main=paste(names(gmtlist)[i], colnames(pg)[j]),
                      Colv=F, Rowv=F, dendrogram="none",
                      colsep=seq(1:dim(plotpg)[[2]]), rowsep=seq(1:dim(plotpg)[[1]]), sepcolor="grey92", sepwidth=c(0.005,0.005),
                      keysize=0.8, col=hcols, density.info="none", symkey=FALSE,
                      labRow=rownames(plotpg), labCol=colnames(plotpg), trace="none",cex.main=0.5, cexRow=cRow, cexCol=0.7, mar=c(10,30) )
            }
        graphics.off()
        pgList[[i]] <- pg; names(pgList)[i] <- names(gmtlist)[i]
     }
return(pgList)
}
###################################################



###################################################
extraGSEAheatmaps <- function(Files, RN, howMany, Color, selectRowsPattern, selPatName, orderBy, print_hmat_xls, histos, pgseamat, onlyInName, matAnn) {
      
      wd <- getwd()
      if( file.exists( "extraGSEAHeatmaps" ) ) { setwd( wd )
        } else { system( paste( "mkdir",  "extraGSEAHeatmaps" , sep=" " ) )
            setwd( wd )
        }
      #Files <- list.files(pattern="Bothsig.+\\.xls$"); KDsigonly <- "BOTHsig"; Files 
      #Files <-  c(list.files(pattern="^GSEA.+GO.xls"),list.files(pattern="^GSEA.+curatedgenesets .xls")); Files; KDsigonly <- "onlyKDsig"; numTop <- 100
      for(j in 1:length(Files)) {
      setwd(wd)              
          FName <- Files[j]
          plist <- read.table( file=paste(wd, Files[j], sep="/"), header=T, stringsAsFactors=F, quote="", fill=T, sep="\t" ); #plist[1:5,1:5]
          colnames(plist)[1] <- "Name"
          plist$L <- as.vector(sapply(plist$genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
      
      setwd( "extraGSEAHeatmaps" )               
          whna <- which( apply(plist[,howMany],1, function(x) any(is.na(x))) )
          if(length(whna)>0) { plist <- plist[-whna ,] }
          
          colNum <- orderBy+1; colN <- colnames(plist)[colNum]
          plist <- plist[order(plist[,colNum]),]
          
          if(selectRowsPattern != "no") { 
                     test <- paste(plist$Name, plist$Description1,  plist$Description2,  plist$Description3,  plist$Description4)
                     if(onlyInName=="yes") { test <- plist$Name }
                     #selectRowsPattern ="BRCA1_OVEREXP_DN|XHM-ED B cells"
                     whpat <- grep(selectRowsPattern, test, ignore.case=T)
                     plist <- plist[whpat,]
                     write.table(plist, file=paste(FName,selPatName, "xls",sep="."), ,row.names=F,quote=F,sep="\t")
          }
          #rn <- dim(plist)[[1]]
          if(RN == "all") { RN <- nrow(plist) }
          rn <- RN
          if(nrow(plist) < rn) { rn <- nrow(plist) }
          plist$L[1:50]
          
          #plot(plist$results.logFCCONvsKD, plist$L)
          #plot(plist$ttest, plist$L)
          #plot(plist$ttest, plist$results.logFCCONvsKD, col="white")
          #text(plist$ttest, plist$results.logFCCONvsKD, labels=plist$L)
          
          #plot(plist$ttest-plist$results.logFCCONvsKD, plist$L)    
          #scatterplot3d(x=plist$ttest, y=plist$results.logFCCONvsKD, z=plist$L, highlight.3d=T)
          
      
      #### genes heatmap ##################################
          pdf(paste("GeneHeatmaps",FName,"pgseamat",rn,"mostSign",selPatName, colN,Color,"pdf", sep="."), height=9, width=12)
          layout(matrix(c(1:12),2,6,byrow=F), widths=rep(2,6), heights=c(9,1), TRUE); #layout.show(12)
          
          if(rn < nrow(plist)/2) { vonWonachWo <- c(1:rn,(nrow(plist)-rn+1):(nrow(plist)) ) }
          if(rn >= nrow(plist)/2) { vonWonachWo <- c(1:rn) }
                                                
          for(i in vonWonachWo ) {
          
                genes <-  as.vector(unlist(strsplit(plist$genes[i],","))); Name <- plist$Name[i]
                whg <- which(rownames(pgseamat) %in% genes)
                pgseamat[whg,]
                #par(mfrow=c(1,2))
          
                #dotchart(pgseamat[whg,1], main=Name, labels=rownames(pgseamat[whg,]),cex.lab=0.5, col=rainbow(length(whg)), pch=20 )
                #abline(v=0, col="grey")
                #dotchart(pgseamat[whg,2], labels=rownames(pgseamat[whg,]),cex.lab=0.5, col=rainbow(length(whg)), pch=20 )
                #abline(v=0, col="grey")
                
                hmat <-  pgseamat[whg,howMany]
                hmat <- hmat[order(hmat[,orderBy]),]
                
                if(histos == "yes") {
                      pdf(file=paste("histos",FName, plist$Name[i], selPatName,"pdf", sep="."), height=9, width=9)
                        par(mfrow=c(2,2))
                           for (xy in 1:length(howMany)) {
                                        plot(density(hmat[,xy]), col="red", main=colnames(hmat)[xy], ylim=c(0,2))
                                        lines(density(pgseamat[,xy]), col="grey")
                                        abline(v=0, col="grey")                          
                                 }
                      dev.off()
                }
                
                if(print_hmat_xls == "yes") {                 
                           whHmat <- which(matAnn$syms %in% rownames(hmat) )
                           write.table(matAnn[whHmat,], file=paste("GeneTable", plist$Name[i], "xls", sep="."), quote=F, row.names=F, sep="\t" )                            
                } 
                
                if(sum(sign(range(hmat, na.rm=T)), na.rm=T) > 0) { hcols <- colorpanel(100, "white", "tomato") } 
                if(sum(sign(range(hmat, na.rm=T)), na.rm=T) < 0) { hcols <- colorpanel(100, "yellowgreen","white") }           
                if(sum(sign(range(hmat, na.rm=T)), na.rm=T) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(hmat),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(hmat), na.rm=T) )) ) }

                if(Color =="BW") {
                    if(sum(sign(range(hmat, na.rm=T)), na.rm=T) > 0) { hcols <- colorpanel(100, "white", "black") } 
                    if(sum(sign(range(hmat, na.rm=T)), na.rm=T) < 0) { hcols <- colorpanel(100, "black","white") }           
                }
          
                #par(cex.main=0.5)
                #heatmap.2(as.matrix(hmat), main=paste(Name),
                #                      Colv=F, Rowv=F, dendrogram="none",
                #                      colsep=seq(1:dim(hmat)[[2]]), rowsep=seq(1:dim(hmat)[[1]]), sepcolor="grey92", sepwidth=c(0.005,0.005),
                #                      keysize=0.8, col=hcols, density.info="none", symkey=FALSE,
                #                      labRow=rownames(hmat), labCol=colnames(hmat), trace="none",cex.main=0.5, cexRow=0.7, cexCol=1, mar=c(10,30) )
                #
                par(mar=c(4, 1, 4, 4)); par(cex.main=0.6)
                m1 <- substr(Name,1,20); m2 <- substr(Name,21,nchar(Name)); 

                image(t(hmat)[,nrow(hmat):1], x=1:ncol(hmat), y=1:nrow(hmat), col=hcols, axes=F, xlab="", ylab="", main=c(m1,m2) )
                axis(1, 1:ncol(hmat), labels=colnames(hmat), cex.axis=0.9, font=2, lwd.ticks=1,las=3)
                axis(4, 1:nrow(hmat), labels=rev(rownames(hmat)), las=2, cex.axis=0.5, lwd.ticks=0)
                abline(v=((1:ncol(hmat))+0.5), col="grey")
                abline(h=((1:nrow(hmat))+0.5), col="grey")
                box(lty = 'solid', col = 'grey3')
                      
                ## draw colorscale on bottom
                par(mar=c(2.5, 2, 4, 4)) 
                matlab <- pretty(sort(as.vector(as.matrix(hmat))),n=10)
                image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                box(lty = 'solid', col = 'black')            
                axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)                                     
                #text(1,1,labels="2", col="black"); text(length(matlab)-0.1,1,labels="13", col="white"); 
          }
          graphics.off()
      }
      setwd(wd)
}
#########################################



#### cluster categories #################
clusterFctCat <- function(Files, whichcols, PvalornumRow, Pval, howMany) {
      
      #whichcols <- c(2:5); PvalornumRow <- "numRow"; Pval <- 0.05; howMany <- 100
      
      for(j in 1:length(Files)) {
          FName <- Files[j]
          tab <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t"); tab[1:5,1:5]
          tab$Name <- gsub(" .+", "", tab$Name)
          tab$L <- as.vector(sapply(tab$genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(tab,1, function(x) any(is.na(x))))
          if(length(whna)>0) { tab <- tab[-whna ,] }
          
          
          if(PvalornumRow == "numRow") {
 
                pdf(file=paste("heatmap",FName,howMany,"pdf",sep="."), height=12, width=9);          
                layout(matrix(c(1:2,3,3,3,3),2,2,byrow=T), widths=c(2,7), heights=c(12,1), TRUE); layout.show(3)                          
                for(wc in 1:length(whichcols)) {
                       tab <- tab[order(tab[,whichcols[wc]]),]
                            
                       ## dmats ########
                       tabn <- tab[1:howMany,]
                       tabp <- tab[(nrow(tab)-howMany+1):nrow(tab),]
                       
                       tabpnList <- list(tabn, tabp)
                       
                       dmatneg <- dmatFct(tabn, "neg")
                       dmatpos <- dmatFct(tabp, "pos")
                       ##################
                                        
                       ## heatmaps ######
                       hcneg <- hclust(as.dist(1-dmatneg),method = "average")
                       hcpos <- hclust(as.dist(1-dmatpos),method = "average")
                       hcList <- list(hcneg, hcpos)
                                              
                       #pdf(file=paste("tmp.neg",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50));  par(cex=0.5);plot(as.dendrogram(hcneg), horiz=T, axes=F, leaflab="perpendicular"); graphics.off()
                       #pdf(file=paste("tmp.pos",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50)); par(cex=0.5); plot(as.dendrogram(hcpos), horiz=T, axes=F); graphics.off()
                              
                       whn <- which(tabn$Name %in% hcneg$labels[rev(hcneg$order)] )
                       whp <- which(tabp$Name %in% hcpos$labels[rev(hcpos$order)] )
                       
                       stopifnot(identical(hcneg$labels[rev(hcneg$order)], tabn$Name[whn][order(tabn$Name)][rev(hcneg$order)] )) 
                       stopifnot(identical(hcpos$labels[rev(hcpos$order)], tabp$Name[whn][order(tabp$Name)][rev(hcpos$order)] )) 
                       
                       for(xx in 1:length(tabpnList)) {                          
                              tabpn <- tabpnList[[xx]]
                              hc <- hcList[[xx]]
                              carp <- t(as.matrix(tabpn[whn[order(tabpn$Name)][rev(hc$order)] ,whichcols]))
                              rn <-   rev(as.character(tabpn[whn[order(tabpn$Name)][rev(hc$order)] ,"Name"]))
                              
                              if(sum(sign(range(carp))) > 0)  { hcols <- colorpanel(100, "white", "tomato") } 
                              if(sum(sign(range(carp))) < 0)  { hcols <- colorpanel(100, "yellowgreen","white") }           
                              if(sum(sign(range(carp))) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(carp),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(carp)))) ) }
                                                                                           
                              par(mar=c(0,0,0,0));plot(as.dendrogram(hcList[[xx]]), horiz=T, axes=F, leaflab="none")                                               
          
                              par(mar=c(2.5,0.1,2.5,30))
                              image(carp, y=1:ncol(carp), x=1:nrow(carp), col=hcols, axes=F, xlab="", ylab="")
                              axis(1, 1:nrow(carp), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                              axis(4, 1:ncol(carp), labels=rn, las=2, cex.axis=0.6, lwd.ticks=0)                          
                              abline(h=((1:ncol(carp))+0.5), col="grey", lwd=0.1)
                              abline(v=((1:nrow(carp))+0.5), col="grey", lwd=0.1)
                              box(lty = 'solid', col = 'grey3')
                       
                             ## draw colorscale on bottom
                             par(mar=c(5, 23, 0, 19)) 
                             matlab <- pretty(sort(as.vector(as.matrix(carp))),n=10)
                             image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                             box(lty = 'solid', col = 'black')            
                             axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1) 
                             #plot(1,1,col="white", axes=F)                       
                       }#; graphics.off()                            
              }                                                      
          }
          graphics.off()
     }
}          
#########################################




#### cluster categories #################
clusterFctCat.v2 <- function(Files, whichcols, whichPcols, PvalornumRow, Pval, howMany) {

      #whichcols <- c(2:5); whichPcols <- c(6:9); PvalornumRow <- "numRow"; Pval <- 0.05; howMany <- 100

      for(j in 1:length(Files)) {
          FName <- Files[j]; print(FName)
          tab <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t"); tab[1:5,1:5]
          tab$Name <- gsub(" .+", "", tab$Name)
          tab$L <- as.vector(sapply(tab$genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(tab,1, function(x) any(is.na(x))))
          if(length(whna)>0) { tab <- tab[-whna ,] }


          if(PvalornumRow == "numRow") {

                pdf(file=paste("heatmap",FName,howMany,"pdf",sep="."), height=12, width=9);
                layout(matrix(c(1,2,3,4,4,5),2,3,byrow=T), widths=c(2,7,2), heights=c(12,1), TRUE); #layout.show(5)
                for(wc in 1:length(whichcols)) {
                       tab <- tab[order(tab[,whichcols[wc]]),]
                       print(paste("wc: ", wc))
                       ## dmats ########
                       tabn <- tab[1:howMany,]
                       tabp <- tab[(nrow(tab)-howMany+1):nrow(tab),]

                       tabpnList <- list(tabn, tabp)

                       dmatneg <- dmatFct(tabn, "neg")
                       dmatpos <- dmatFct(tabp, "pos")
                       ##################

                       ## heatmaps ######
                       hcneg <- hclust(as.dist(1-dmatneg),method = "average")
                       hcpos <- hclust(as.dist(1-dmatpos),method = "average")
                       hcList <- list(hcneg, hcpos)

                       #pdf(file=paste("tmp.neg",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50));  par(cex=0.5);plot(as.dendrogram(hcneg), horiz=T, axes=F, leaflab="perpendicular"); graphics.off()
                       #pdf(file=paste("tmp.pos",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50)); par(cex=0.5); plot(as.dendrogram(hcpos), horiz=T, axes=F); graphics.off()

                       whn <- which(tabn$Name %in% hcneg$labels[rev(hcneg$order)] )
                       whp <- which(tabp$Name %in% hcpos$labels[rev(hcpos$order)] )

                       stopifnot(identical(hcneg$labels[rev(hcneg$order)], tabn$Name[whn][order(tabn$Name)][rev(hcneg$order)] ))
                       stopifnot(identical(hcpos$labels[rev(hcpos$order)], tabp$Name[whp][order(tabp$Name)][rev(hcpos$order)] ))

#pdf(file=paste("heatmap",FName,howMany,"pdf",sep="."), height=12, width=9);
#layout(matrix(c(1,2,3,4,4,5),2,3,byrow=T), widths=c(2,7,2), heights=c(12,1), TRUE); #layout.show(5)
                       for(xx in 1:length(tabpnList)) {
                              tabpn <- tabpnList[[xx]]
                              hc <- hcList[[xx]]
                              carp <-  t(as.matrix(tabpn[whn[order(tabpn$Name)][rev(hc$order)] ,whichcols]))
                              carpP <- t(as.matrix(tabpn[whn[order(tabpn$Name)][rev(hc$order)] ,whichPcols])); carpP <- -log(carpP,10)
                              rn <-   rev(as.character(tabpn[whn[order(tabpn$Name)][rev(hc$order)] ,"Name"]))

                              ## carpet of z-Values ################
                              if(sum(sign(range(carp))) > 0) { hcols <- colorpanel(100, "white", "tomato") }
                              if(sum(sign(range(carp))) < 0) { hcols <- colorpanel(100, "yellowgreen","white") }
                              if(sum(sign(range(carp))) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(carp),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(carp)))) ) }

                              par(mar=c(0,0,0,0));plot(as.dendrogram(hcList[[xx]]), horiz=T, axes=F, leaflab="none")

                              par(mar=c(2.6, 0.1, 2.6, 35))
                              image(carp, y=1:ncol(carp), x=1:nrow(carp), col=hcols, axes=F, xlab="", ylab="")
                              axis(1, 1:nrow(carp), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                              axis(4, 1:ncol(carp), labels=rn, las=2, cex.axis=0.6, lwd.ticks=0)
                              abline(h=((1:ncol(carp))+0.5), col="grey", lwd=0.1)
                              abline(v=((1:nrow(carp))+0.5), col="grey", lwd=0.1)
                              box(lty = 'solid', col = 'grey3')
                              ######################################
                              
                              ## carpet of P values ################
                              Pcols <- colorpanel(100, "honeydew", "yellow","magenta2")
                              par(mar=c(2.6, 4, 2.6, 3))
                              image(carpP, y=1:ncol(carpP), x=1:nrow(carpP), col=Pcols, axes=F, xlab="", ylab="")
                              axis(1, 1:nrow(carpP), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                              #axis(4, 1:ncol(carpP), labels="", las=2, cex.axis=0.6, lwd.ticks=0)
                              abline(h=((1:ncol(carpP))+0.5), col="grey", lwd=0.1)
                              abline(v=((1:nrow(carpP))+0.5), col="grey", lwd=0.1)
                              box(lty = 'solid', col = 'grey3')
                              #######################################

                              ## draw z-Value colorscale on bottom ##
                              par(mar=c(2, 10, 3.5, 35))
                              matlab <- pretty(sort(as.vector(as.matrix(carp))),n=10)
                              image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                              box(lty = 'solid', col = 'black')
                              axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)
                              #######################################
                              
                              ## draw P-Value colorscale on bottom ##
                              par(mar=c(2, 2, 3.5, 3))
                              matlab <- pretty(sort(as.vector(as.matrix(carpP))),n=6)
                              image(as.matrix(matlab), x=1:length(matlab),y=1, col=Pcols, axes=F, xlab="", ylab="")
                              box(lty = 'solid', col = 'black')
                              axis(1, 1:length(matlab), labels=as.character(round(10^-matlab,3)), las=3, cex.axis=0.5, lwd.ticks=1)
                              #######################################
                              
                       }; #graphics.off()
              }
          }
          graphics.off()
     }
}
#########################################





### dmat Fct ############################
dmatFct <- function(Tab, Name) {
     dmat <- matrix(ncol=length(levels(as.factor(Tab$Name))), nrow=length(levels(as.factor(Tab$Name)))); dmat              
     rownames(dmat) <- colnames(dmat) <- levels(as.factor(Tab$Name))
     dmatgenes <- as.data.frame(dmat,stringsAsFactors=F)
  
     for (l in 1:length(levels(as.factor(Tab$Name))) ) {
          for (m in 1:length(levels(as.factor(Tab$Name))) ) {
               wh1 <- levels(as.factor(Tab$Name))[l]
               wh2 <- levels(as.factor(Tab$Name))[m]
               a <- as.vector(unlist(strsplit(Tab[Tab$Name == wh1,]$genes,","))); a <- a[which(!is.na(a))]
               b <- as.vector(unlist(strsplit(Tab[Tab$Name == wh2,]$genes,","))); b <- b[which(!is.na(b))]
               c <- intersect(a,b); c
               myDist <- length(c)/(min(length(a),length(b))); myDist
               dmat[l,m] <- myDist; #colnames(dmat)[i] <- levels(as.factor(Tab$Name))[i]
               dmatgenes[l,m] <- paste(c,collapse=",")
           }
     }
     dmat[dmat == 0] <- 0.00001
     dmat[dmat == 1] <- 0.99999                 
assign(paste("dmat",Name,sep=""),dmat)
Func.env$writeOutxls(dmatgenes,Name)
dmat
}
#########################################



### dmat Fct DAVID ######################
dmatFctDavid <- function(Tab, Name) {
     dmat <- matrix(ncol=length(levels(as.factor(Tab$Term))), nrow=length(levels(as.factor(Tab$Term)))); dmat              
     rownames(dmat) <- colnames(dmat) <- levels(as.factor(Tab$Term))
  
     for (l in 1:length(levels(as.factor(Tab$Term))) ) {
          for (m in 1:length(levels(as.factor(Tab$Term))) ) {
               wh1 <- levels(as.factor(Tab$Term))[l]
               wh2 <- levels(as.factor(Tab$Term))[m]
               a <- as.vector(unlist(strsplit(Tab[Tab$Term == wh1,]$Genes,","))); a <- a[which(!is.na(a))]
               b <- as.vector(unlist(strsplit(Tab[Tab$Term == wh2,]$Genes,","))); b <- b[which(!is.na(b))]
               dmat[l,m] <- length(intersect(a,b)) / length(union(a,b))
           }
     }
     dmat[dmat == 0] <- 0.00001
     dmat[dmat == 1] <- 0.99999                 
assign(paste("dmat",Name,sep=""),dmat)
}
#########################################




#### cluster categories #################
clusterFctCat.v3 <- function(Files, whichcols, whichPcols, PvalornumRow, Pval, howMany) {

      #whichcols <- c(2:5); whichPcols <- c(6:9); PvalornumRow <- "Pval"; Pval <- 0.05; howMany <- 100

      for(j in 1:length(Files)) {
          FName <- Files[j]; print(FName)
          tab <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t"); tab[1:5,1:5]
          tab$Name <- gsub(" .+", "", tab$Name)
          tab$L <- as.vector(sapply(tab$genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(tab,1, function(x) any(is.na(x))))
          if(length(whna)>0) { tab <- tab[-whna ,] }

          if(PvalornumRow == "numRow") {
                 pdf(file=paste("heatmap",FName,"numRow",howMany,"pdf",sep="."), height=12, width=9);
          }
          if(PvalornumRow == "Pval") {
                 pdf(file=paste("heatmap",FName,"Pval",Pval,"pdf",sep="."), height=12, width=9);
          }
          layout(matrix(c(1,2,3,4,4,5),2,3,byrow=T), widths=c(2,7,2), heights=c(12,1), TRUE); #layout.show(5)


          for(wc in 1:length(whichcols)) {
                 if(PvalornumRow == "numRow") {
                       tab <- tab[order(tab[,whichcols[wc]]),]
                       print(paste("wc: ", wc))
                       tabn <- tab[1:howMany,]
                       tabp <- tab[(nrow(tab)-howMany+1):nrow(tab),]
                  }
                 if(PvalornumRow == "Pval") {
                       tabn <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] < 0) ), ]
                       tabp <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] > 0) ), ]
                  }
                  
                       tabpnList <- list(tabn, tabp)

                       dmatneg <- dmatFct(tabn, "neg")
                       dmatpos <- dmatFct(tabp, "pos")
                       ##################

                       ## heatmaps ######
                       hcneg <- hclust(as.dist(1-dmatneg),method = "average")
                       hcpos <- hclust(as.dist(1-dmatpos),method = "average")
                       hcList <- list(hcneg, hcpos)

                       #pdf(file=paste("tmp.neg",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50));  par(cex=0.5);plot(as.dendrogram(hcneg), horiz=T, axes=F, leaflab="perpendicular"); graphics.off()
                       #pdf(file=paste("tmp.pos",howMany,"col",wc,"pdf",sep="."), height=12, width=9); par(mar=c(2,0,0,50)); par(cex=0.5); plot(as.dendrogram(hcpos), horiz=T, axes=F); graphics.off()

                       whn <- which(tabn$Name %in% hcneg$labels[rev(hcneg$order)] )
                       whp <- which(tabp$Name %in% hcpos$labels[rev(hcpos$order)] )

                       stopifnot(identical(hcneg$labels[rev(hcneg$order)], tabn$Name[whn][order(tabn$Name)][rev(hcneg$order)] ))
                       stopifnot(identical(hcpos$labels[rev(hcpos$order)], tabp$Name[whp][order(tabp$Name)][rev(hcpos$order)] ))

#pdf(file=paste("heatmap",FName,howMany,"pdf",sep="."), height=12, width=9);
#layout(matrix(c(1,2,3,4,4,5),2,3,byrow=T), widths=c(2,7,2), heights=c(12,1), TRUE); #layout.show(5)
                       for(xx in 1:length(tabpnList)) {
                              tabpn <- tabpnList[[xx]]
                              hc <- hcList[[xx]]
                              wh <- which(tabpn$Name %in% hc$labels[rev(hc$order)] )

                              carp <-  t(as.matrix(tabpn[wh[order(tabpn$Name)][rev(hc$order)] ,whichcols]))
                              carpP <- t(as.matrix(tabpn[wh[order(tabpn$Name)][rev(hc$order)] ,whichPcols])); carpP <- -log(carpP,10)
                              rn <-   rev(as.character(tabpn[wh[order(tabpn$Name)][rev(hc$order)] ,"Name"]))

                              ## carpet of z-Values ################
                              if(sum(sign(range(carp))) > 0) { hcols <- colorpanel(100, "white", "tomato") }
                              if(sum(sign(range(carp))) < 0) { hcols <- colorpanel(100, "yellowgreen","white") }
                              if(sum(sign(range(carp))) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(carp),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(carp)))) ) }

                              par(mar=c(0,0,0,0));plot(as.dendrogram(hcList[[xx]]), horiz=T, axes=F, leaflab="none")

                              par(mar=c(2.6, 0.1, 2.6, 35))
                              image(carp, y=1:ncol(carp), x=1:nrow(carp), col=hcols, axes=F, xlab="", ylab="")
                              axis(1, 1:nrow(carp), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                              axis(4, 1:ncol(carp), labels=rn, las=2, cex.axis=0.6, lwd.ticks=0)
                              abline(h=((1:ncol(carp))+0.5), col="grey", lwd=0.1)
                              abline(v=((1:nrow(carp))+0.5), col="grey", lwd=0.1)
                              box(lty = 'solid', col = 'grey3')
                              ######################################

                              ## carpet of P values ################
                              Pcols <- colorpanel(100, "honeydew", "yellow","magenta2")
                              par(mar=c(2.6, 4, 2.6, 3))
                              image(carpP, y=1:ncol(carpP), x=1:nrow(carpP), col=Pcols, axes=F, xlab="", ylab="")
                              axis(1, 1:nrow(carpP), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                              #axis(4, 1:ncol(carpP), labels="", las=2, cex.axis=0.6, lwd.ticks=0)
                              abline(h=((1:ncol(carpP))+0.5), col="grey", lwd=0.1)
                              abline(v=((1:nrow(carpP))+0.5), col="grey", lwd=0.1)
                              box(lty = 'solid', col = 'grey3')
                              #######################################

                              ## draw z-Value colorscale on bottom ##
                              par(mar=c(2, 10, 3.5, 35))
                              matlab <- pretty(sort(as.vector(as.matrix(carp))),n=10)
                              image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                              box(lty = 'solid', col = 'black')
                              axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)
                              #######################################

                              ## draw P-Value colorscale on bottom ##
                              par(mar=c(2, 2, 3.5, 3))
                              matlab <- pretty(sort(as.vector(as.matrix(carpP))),n=6)
                              image(as.matrix(matlab), x=1:length(matlab),y=1, col=Pcols, axes=F, xlab="", ylab="")
                              box(lty = 'solid', col = 'black')
                              axis(1, 1:length(matlab), labels=as.character(round(10^-matlab,3)), las=3, cex.axis=0.5, lwd.ticks=1)
                              #######################################

                       }; #graphics.off()
              }
          graphics.off()
     }
}
#########################################





