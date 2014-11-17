#### cluster categories #################
clusterFctCatDAVID.v3 <- function(Files, whichcols, whichPcols, PvalornumRow, Pval, howmany, largeCat, Plot, carps, pdftif, cut.col.scheme, col.overRide="") {

   carpList <- list() 
   for(j in 1:length(Files)) {
          
          #whichcols <- c(12); whichPcols <- c(12); PvalornumRow <- "numRow"; Pval <- 0.05; howmany <- numCat; largeCat <- 1000; Plot <- "NO"; carps=NULL; pdftif=PT; cut.col.scheme=0.4; col.overRide=""
           
          FName <- Files[j]; print(FName)
          tab <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t", quote = "\""); tab[1:5,1:5]
		  tab$Category <- gsub("GOTERM_BP_FAT", "GO", tab$Category)
		  tab$Category <- gsub("GOTERM_MF_FAT", "GO", tab$Category)
		  tab$Category <- gsub("GOTERM_CC_FAT", "GO", tab$Category)
		  tab$Term <- gsub("GO:", "", tab$Term)
		  tab$Term <- paste(tolower(tab$Category), tab$Term, sep=":")
          howMany <- min(howmany,dim(tab)[[1]])


          tab$L <- as.vector(sapply(tab$Genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(tab,1, function(x) any(is.na(x))))
          if(length(whna)>0) { tab <- tab[-whna ,] }

          whL <- which(tab$Pop.Hits >= largeCat)
          if(length(whL)>0) { tab <- tab[-whL,] }
          
          if(PvalornumRow == "numRow") {
                if(pdftif == "pdf") {
                    pdf(file=paste(FName,"heatmap","numRow",howMany,"largeCatexcl",largeCat,"pdf","part",sep="."), height=12, width=9);
                }
                if(pdftif == "tiff") {             
                    tiff(filename=paste(FName,"heatmap","numRow",howMany,"largeCatexcl",largeCat,"tiff","part",sep="."), height=12, width=9,
                                 units="in", pointsize=12, compression="lzw", bg="white", res=700)           
                }
           
           }
           if(PvalornumRow == "Pval") {
             pdf(file=paste(FName,"heatmap","Pval",Pval,"largeCatexcl",largeCat,"col.scheme",col.scheme,"pdf","part",sep="."), height=12, width=9);
           }
          
          layout(matrix(c(1:2,3,3,3,3),2,2,byrow=T), widths=c(2,7), heights=c(12,1), TRUE); #layout.show(3)                          


          for(wc in 1:length(whichcols)) {
               if(PvalornumRow == "numRow") {
                     tab <- tab[order(tab[,whichcols[wc]]),]
                     print(paste("wc: ", wc))
                     tabn <- tab[1:howMany,]
                     #tabp <- tab[(nrow(tab)-howMany+1):nrow(tab),]
                
                     numlow <- length(which(tabn[,whichPcols] > 0.3))
                     if(numlow >= howMany*cut.col.scheme) { col.scheme <- "lowsign" }
                     if(numlow < howMany*cut.col.scheme) { col.scheme <- "highsign" }
                     if(nchar(col.overRide) > 1) { col.scheme <- col.overRide }
                                     
                }
               if(PvalornumRow == "Pval") {
                     tabn <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] < 0) ), ]
                     #tabp <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] > 0) ), ]
                }
                

               dmat <- dmatFctDavid(tabn, "top")
               hc <- hclust(as.dist(1-dmat),method = "average")
               #pdf(file="testdendro.pdf", height=12, width=9)
               #      par(mar=c(0,0,0,30))
               #      plot(as.dendrogram(hc), horiz=T)
               #graphics.off()
               
               whn <- which(tabn$Term %in% hc$labels[rev(hc$order)] )
               stopifnot(identical(hc$labels[rev(hc$order)], tabn$Term[whn][order(tabn$Term)][rev(hc$order)] )) # if this assertion fails, check if you have duplicated term names in your list

                wh <- which(tabn$Term %in% hc$labels[rev(hc$order)] )

                #carp <-  t(as.matrix(tabn[wh[order(tabn$Term)][rev(hc$order)] ,whichcols]))
                carpP <- t(as.matrix(rev(tabn[wh[order(tabn$Term)][rev(hc$order)] ,whichPcols]))); carpP <- -log(carpP,10)
                rn <-   rev(as.character(tabn[wh[order(tabn$Term)][rev(hc$order)] ,"Term"]))

                #if (j==1) { CARP <- carpP }
                #if (j>1) { CARP <- rbind(CARP,carpP) }
                carpList[[j]] <- carpP; print(c(j, carpP))
                if( Plot == "NO" ) { next; }

                if( Plot == "YES" ) {
                
                      ## carpet of z-Values ################
                      #if(sum(sign(range(carp))) > 0) { hcols <- colorpanel(100, "white", "tomato") }
                      #if(sum(sign(range(carp))) < 0) { hcols <- colorpanel(100, "yellowgreen","white") }
                      #if(sum(sign(range(carp))) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(carp),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(carp)))) ) }
      
                      par(mar=c(0,0,0,0));plot(as.dendrogram(hc), horiz=T, axes=F, leaflab="none")
      
                      #par(mar=c(2.6, 0.1, 2.6, 35))
                      #image(carp, y=1:ncol(carp), x=1:nrow(carp), col=hcols, axes=F, xlab="", ylab="")
                      #axis(1, 1:nrow(carp), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                      #axis(4, 1:ncol(carp), labels=rn, las=2, cex.axis=0.6, lwd.ticks=0)
                      #abline(h=((1:ncol(carp))+0.5), col="grey", lwd=0.1)
                      #abline(v=((1:nrow(carp))+0.5), col="grey", lwd=0.1)
                      #box(lty = 'solid', col = 'grey3')
                      ######################################
      
                      ## carpet of P values ################
                      #allcols <- colorpanel(length(carps[[j]]), "honeydew", "green","magenta2")
                      allcols <- sort(carps[[j]]); allcols1 <- 10^-allcols
                      colvec <- c("grey45", "grey55", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red", "VioletRed1", "DeepPink", "DeepPink2")
                      matlab <- c(">0.5", 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001, 0.00005)

                      allcols[which(allcols1 >= 0.5)] <- "grey45"
                      allcols[which((allcols1  < 0.5) & (allcols1 >=0.4))] <- "grey55"
                      allcols[which((allcols1  < 0.4) & (allcols1 >=0.3))] <- "SlateGray3"
                      allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "SlateGray1"
                      allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "DarkSeaGreen1"
                      allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "yellow"
                      allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "orange"
                      allcols[which((allcols1  < 0.01) & (allcols1 >=0.005))] <- "tomato"
                      allcols[which((allcols1  < 0.005) & (allcols1 >=0.001))] <- "red"
                      allcols[which((allcols1  < 0.001) & (allcols1 >=0.0001))] <- "VioletRed1"
                      allcols[which((allcols1  < 0.0001) & (allcols1 >=0.00005))] <- "DeepPink"
                      allcols[which((allcols1  < 0.00005))] <- "DeepPink2"
                      
                      if (col.scheme == "lowsign") {
                            colvec <- c("grey10", "grey25", "grey40", "grey60", "grey80", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red")
                            matlab <- c(">0.9", 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001)
      
                            allcols[which(allcols1 >= 0.8)] <- "grey10"
                            allcols[which((allcols1  < 0.8) & (allcols1 >=0.7))] <- "grey25"
                            allcols[which((allcols1  < 0.7) & (allcols1 >=0.6))] <- "grey40"
                            allcols[which((allcols1  < 0.6) & (allcols1 >=0.5))] <- "grey60"
                            allcols[which((allcols1  < 0.5) & (allcols1 >=0.4))] <- "grey80"
                            allcols[which((allcols1  < 0.4) & (allcols1 >=0.3))] <- "SlateGray3"
                            allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "SlateGray1"
                            allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "DarkSeaGreen1"
                            allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "yellow"
                            allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "orange"
                            allcols[which((allcols1  < 0.01) & (allcols1 >=0.005))] <- "tomato"
                            allcols[which((allcols1  < 0.005) & (allcols1 >=0.001))] <- "VioletRed1"                            
                            allcols[which((allcols1  < 0.001) & (allcols1 >=0.0001))] <- "DeepPink"                                                        
                            allcols[which((allcols1  < 0.0001))] <- "DeepPink2"                                            
                      }
                      
                      if (col.scheme == "highsign") {
                      colvec <- c("grey45", "grey55", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red", "VioletRed1", "DeepPink", "DeepPink2")
                            matlab <- c(">0.3", 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001)
      
                            allcols[which(allcols1 >= 0.3)] <- "grey45"
                            allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "grey55"
                            allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "SlateGray3"
                            allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "SlateGray1"
                            allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "DarkSeaGreen1"
                            allcols[which((allcols1  < 0.01) & (allcols1 >=0.001))] <- "yellow"
                            allcols[which((allcols1  < 0.001) & (allcols1 >=0.0001))] <- "orange"
                            allcols[which((allcols1  < 0.0001) & (allcols1 >=0.00001))] <- "tomato"
                            allcols[which((allcols1  < 0.00001) & (allcols1 >=0.000001))] <- "red"
                            allcols[which((allcols1  < 0.000001) & (allcols1 >=0.0000001))] <- "VioletRed1"
                            allcols[which((allcols1  < 0.0000001) & (allcols1 >=0.00000001))] <- "DeepPink"
                            allcols[which((allcols1  < 0.000000001))] <- "DeepPink2"                                            
                      }                      
                      
                      
                      #pie(rep(1,length(allcols)), col=allcols)                      
                      
                      whcol <- which(as.character(sort(carps[[j]])) %in% as.character(carpP))
                      Pcols <- allcols[whcol]
                      #pie(rep(1,length(Pcols)), col=Pcols)
                      
                      ## um dopplete zu vermeiden ##
                      or <- order(carpP,decreasing=F)
                      r <-  vector(length=length(carpP))
                      r[or] <- c(1:length(carpP))
                      ##############################

                      #######################################
                      par(mar=c(2.5,0.1,2.5,30))
                      image(matrix(r,nrow=1), y=1:ncol(carpP), x=1:nrow(carpP), col=Pcols, axes=F, xlab="", ylab="", main=FName, cex.main=0.6)
                      axis(1, 1:nrow(carpP), labels="", las=2, cex.axis=1.5, lwd.ticks=0)
                      axis(4, 1:ncol(carpP), labels=rn, las=2, cex.axis=0.8, lwd.ticks=0)
                      abline(h=((1:ncol(carpP))+0.5), col="grey", lwd=0.1)
                      abline(v=((1:nrow(carpP))+0.5), col="grey", lwd=0.1)
                      box(lty = 'solid', col = 'grey3')
                      #######################################
      
                      ## draw z-Value colorscale on bottom ##
                      #par(mar=c(2, 10, 3.5, 35))
                      #matlab <- pretty(sort(as.vector(as.matrix(carp))),n=10)
                      #image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                      #box(lty = 'solid', col = 'black')
                      #axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)
                      #######################################
      
                      ## draw P-Value colorscale on bottom ##
                      par(mar=c(4, 0, 0, 25))  #par(mar=c(2, 2, 3.5, 3))
                      image(as.matrix(c(1:length(matlab))), x=1:length(matlab),y=1, col=colvec, axes=F, xlab="", ylab="")
                      box(lty = 'solid', col = 'black')
                      axis(1, 1:length(matlab), labels=as.character(matlab), las=3, cex.axis=0.8, lwd.ticks=1)
                       
                      #matlab <- pretty(sort(as.vector(as.matrix(carps[[j]]))),n=6)
                      #image(as.matrix(matlab), x=1:length(matlab),y=1, col=allcols, axes=F, xlab="", ylab="")
                      #box(lty = 'solid', col = 'black')
                      #axis(1, 1:length(matlab), labels=as.character(round(10^-matlab,3)), las=3, cex.axis=0.5, lwd.ticks=1)
                      #######################################
                      #graphics.off()
                 }
          }
      graphics.off()
     }
  if( Plot == "NO" ) { carpList }
}
###################################################





###################################################
DAVIDGeneHeatmaps <- function(Files, RN, howMany, pmat, Color) {

      #Files <- list.files(pattern="Bothsig.+\\.xls$"); KDsigonly <- "BOTHsig"; Files
      #Files <-  c(list.files(pattern="^GSEA.+GO.xls"),list.files(pattern="^GSEA.+curatedgenesets .xls")); Files; KDsigonly <- "onlyKDsig"; numTop <- 100
      for(j in 1:length(Files)) {
          FName <- Files[j]
          plist <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t", quote = "\""); plist[1:5,1:5]
          plist$L <- as.vector(sapply(plist$Genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(plist,1, function(x) any(is.na(x))))
          if(length(whna)>0) { plist <- plist[-whna ,] }

          colNum <- 12; colN <- colnames(plist)[colNum]
          plist <- plist[order(plist[,colNum]),]

          #rn <- dim(plist)[[1]]
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
          pdf(paste("Geneheatmap",FName,"DAVID",rn,"mostSignCategories",colN,"Color",Color,"pdf", sep="."), height=9, width=12)
          layout(matrix(c(1:12),2,6,byrow=F), widths=rep(2,6), heights=c(9,1), TRUE); #layout.show(12)
          for(i in c(1:rn,(nrow(plist)-rn+1):(nrow(plist)) ) ) {
                genes <-  as.vector(unlist(strsplit(plist$Genes[i],","))); Name <- plist$Term[i]
                genes <- gsub(" ", "", genes)
                whg <- which(rownames(pmat) %in% genes)
                pmat[whg,]
                #par(mfrow=c(1,2))

                #dotchart(pmat[whg,1], main=Name, labels=rownames(pmat[whg,]),cex.lab=0.5, col=rainbow(length(whg)), pch=20 )
                #abline(v=0, col="grey")
                #dotchart(pmat[whg,2], labels=rownames(pmat[whg,]),cex.lab=0.5, col=rainbow(length(whg)), pch=20 )
                #abline(v=0, col="grey")

                hmat <-  pmat[whg,howMany]
                hmat <- as.matrix(hmat[order(hmat[,1]),])
                if(dim(hmat)[[1]] == 1) { hmat <- rbind(hmat,hmat) }

                if(sum(sign(range(hmat))) > 0) { hcols <- colorpanel(10, "white", "tomato") }
                if(sum(sign(range(hmat))) < 0) { hcols <- colorpanel(10, "yellowgreen","white") }
                if(sum(sign(range(hmat))) == 0){ hcols1 <- SLmisc$heatmapCol(data=as.matrix(hmat),  col=colorpanel(10, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(hmat)))) ) 
                                                  
                                if(length(hmat[hmat > 0]) > length(hmat[hmat < 0]) ) {                   
                                                  maxNull <- max(which(hcols1 == "#FFFFFF"))
                                                  hcols2 <- colorpanel(length(hcols1)-maxNull+1, "white", "tomato")
                                                  hcols <- hcols1                                                  
                                                  hcols[maxNull:length(hcols)] <- hcols2
                                 }                 
                                if(length(hmat[hmat < 0]) > length(hmat[hmat > 0]) ) {                   
                                                  minNull <- min(which(hcols1 == "#FFFFFF"))
                                                  hcols2 <- colorpanel(minNull, "yellowgreen", "white")
                                                  hcols <- hcols1                                                  
                                                  hcols[1:minNull] <- hcols2
                                 }                                                   
                
                }
                # pdf(); pie(rep(1,length(hcols)), col=hcols,  cex.lab=0.4 ); graphics.off()
                
                if(Color == "BW") {
                    if(sum(sign(range(hmat))) > 0) { hcols <- colorpanel(10, "white", "black") }
                    if(sum(sign(range(hmat))) < 0) { hcols <- colorpanel(10, "black","white") }
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
                #matlab[1] <- min(sort(hmat)); matlab[length(matlab)] <- max(sort(hmat)); 
                #matlab <- round(matlab,1)
                #whh <- which(round(sort(hmat),1) %in% matlab)
                #whdup <- which(duplicated(round(sort(hmat),1)))
                #wh <- setdiff(whh, intersect(whh, whdup))
                                
                image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                box(lty = 'solid', col = 'black')
                axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)
                #text(1,1,labels="2", col="black"); text(length(matlab)-0.1,1,labels="13", col="white");
          }
          graphics.off()
      }
}
#########################################



#### cluster categories #################
clusterFctCatDAVID.specOrd.v2 <- function(Files, whichcols, whichPcols, PvalornumRow, Pval, howmany, largeCat, Plot, carps, pdftif, cut.col.scheme,  newFctCh) {

   carpList <- list() 
   for(j in 1:length(Files)) {
          
          #whichcols <- c(12); whichPcols <- c(12); PvalornumRow <- "numRow"; Pval <- 0.05; howmany <- numCat; largeCat <- 10000; Plot <- "NO"; carps=NULL; pdftif=PT
           
          FName <- Files[j]; print(FName)
          tab <- read.table(file=Files[j], header=T, stringsAsFactors=F, fill=T, sep="\t", quote = "\""); tab[1:5,1:5]
          howMany <- min(howmany,dim(tab)[[1]])


          tab$L <- as.vector(sapply(tab$Genes,function(x) length(unlist(strsplit(x,","))), simplify=T))
          whna <- which(apply(tab,1, function(x) any(is.na(x))))
          if(length(whna)>0) { tab <- tab[-whna ,] }

          whL <- which(tab$Pop.Hits >= largeCat)
          if(length(whL)>0) { tab <- tab[-whL,] }
          
          if(PvalornumRow == "numRow") {
                if(pdftif == "pdf") {
                    pdf(file=paste("heatmap",FName,"numRow",howMany,"largeCatexcl",largeCat,"pdf",sep="."), height=12, width=9);
                }
                if(pdftif == "tiff") {             
                    tiff(filename=paste("heatmap",FName,"numRow",howMany,"largeCatexcl",largeCat,"tiff",sep="."), height=12, width=9,
                                 units="in", pointsize=12, compression="lzw", bg="white", res=700)           
                }
           
           }
           if(PvalornumRow == "Pval") {
             pdf(file=paste("heatmap",FName,"Pval",Pval,"largeCatexcl",largeCat,"col.scheme",col.scheme,"pdf",sep="."), height=12, width=9);
           }
          
          layout(matrix(c(1:2,3,3,3,3),2,2,byrow=T), widths=c(2,7), heights=c(12,1), TRUE); #layout.show(3)                          


          for(wc in 1:length(whichcols)) {
               if(PvalornumRow == "numRow") {
                     tab <- tab[order(tab[,whichcols[wc]]),]
                     print(paste("wc: ", wc))
                     tabn <- tab[1:howMany,]
                     #tabp <- tab[(nrow(tab)-howMany+1):nrow(tab),]
                
                     numlow <- length(which(tabn[,whichPcols] > 0.3))
                     if(numlow >= howMany*cut.col.scheme) { col.scheme <- "lowsign" }
                     if(numlow < howMany*cut.col.scheme) { col.scheme <- "highsign" }                
                }
               if(PvalornumRow == "Pval") {
                     tabn <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] < 0) ), ]
                     #tabp <- tab[intersect(which(tab[whichPcols[wc]] <= Pval), which(tab[whichcols[wc]] > 0) ), ]
                }
                

               dmat <- dmatFctDavid(tabn, "top")
               hc <- hclust(as.dist(1-dmat),method = "average")
               
               ### special part for Gerhards groups ##
               orderNew <- hc$order
               for(xx in 1:length(Lev)) {
                       whLev <- which(newFctCh$Category.KD == Lev[xx])
                       whTerm <-which(hc$labels %in% newFctCh$Term[whLev])
                       origOrder <- which(hc$order %in% whTerm) 
                       Termor <- order(newFctCh$Benjamini[whTerm])
                       
                       whTerm <- whTerm[Termor]
                       for(xxx in 1:length(whTerm)) {
                             orderNew[origOrder[Termor[xxx]] ]    <- whTerm[Termor[xxx]]                       
                       }
               }
               hc$order <- orderNew
               #######################################
               
               
               #pdf(file="testdendro.pdf", height=12, width=9)
               #      par(mar=c(0,0,0,30))
               #      plot(as.dendrogram(hc), horiz=T)
               #graphics.off()
               
               whn <- which(tabn$Term %in% hc$labels[rev(hc$order)] )
               stopifnot(identical(hc$labels[rev(hc$order)], tabn$Term[whn][order(tabn$Term)][rev(hc$order)] ))

                wh <- which(tabn$Term %in% hc$labels[rev(hc$order)] )

                #carp <-  t(as.matrix(tabn[wh[order(tabn$Term)][rev(hc$order)] ,whichcols]))
                carpP <- t(as.matrix(rev(tabn[wh[order(tabn$Term)][rev(hc$order)] ,whichPcols]))); carpP <- -log(carpP,10)
                rn <-   rev(as.character(tabn[wh[order(tabn$Term)][rev(hc$order)] ,"Term"]))

                #if (j==1) { CARP <- carpP }
                #if (j>1) { CARP <- rbind(CARP,carpP) }
                carpList[[j]] <- carpP; print(c(j, carpP))
                if( Plot == "NO" ) { next; }

                if( Plot == "YES" ) {
                
                      ## carpet of z-Values ################
                      #if(sum(sign(range(carp))) > 0) { hcols <- colorpanel(100, "white", "tomato") }
                      #if(sum(sign(range(carp))) < 0) { hcols <- colorpanel(100, "yellowgreen","white") }
                      #if(sum(sign(range(carp))) == 0) { hcols <- SLmisc$heatmapCol(data=as.matrix(carp),  col=colorpanel(100, "yellowgreen","white", "tomato") , lim=min(abs(range(as.matrix(carp)))) ) }
      
                      par(mar=c(0,0,0,0));plot(as.dendrogram(hc), horiz=T, axes=F, leaflab="none")
      
                      #par(mar=c(2.6, 0.1, 2.6, 35))
                      #image(carp, y=1:ncol(carp), x=1:nrow(carp), col=hcols, axes=F, xlab="", ylab="")
                      #axis(1, 1:nrow(carp), labels=rownames(carp), las=2, cex.axis=1.5, lwd.ticks=0)
                      #axis(4, 1:ncol(carp), labels=rn, las=2, cex.axis=0.6, lwd.ticks=0)
                      #abline(h=((1:ncol(carp))+0.5), col="grey", lwd=0.1)
                      #abline(v=((1:nrow(carp))+0.5), col="grey", lwd=0.1)
                      #box(lty = 'solid', col = 'grey3')
                      ######################################
      
                      ## carpet of P values ################
                      #allcols <- colorpanel(length(carps[[j]]), "honeydew", "green","magenta2")
                      allcols <- sort(carps[[j]]); allcols1 <- 10^-allcols
                      colvec <- c("grey45", "grey55", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red", "VioletRed1", "DeepPink", "DeepPink2")
                      matlab <- c(">0.5", 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001, 0.00005)

                      allcols[which(allcols1 >= 0.5)] <- "grey45"
                      allcols[which((allcols1  < 0.5) & (allcols1 >=0.4))] <- "grey55"
                      allcols[which((allcols1  < 0.4) & (allcols1 >=0.3))] <- "SlateGray3"
                      allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "SlateGray1"
                      allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "DarkSeaGreen1"
                      allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "yellow"
                      allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "orange"
                      allcols[which((allcols1  < 0.01) & (allcols1 >=0.005))] <- "tomato"
                      allcols[which((allcols1  < 0.005) & (allcols1 >=0.001))] <- "red"
                      allcols[which((allcols1  < 0.001) & (allcols1 >=0.0001))] <- "VioletRed1"
                      allcols[which((allcols1  < 0.0001) & (allcols1 >=0.00005))] <- "DeepPink"
                      allcols[which((allcols1  < 0.00005))] <- "DeepPink2"
                      
                      if (col.scheme == "lowsign") {
                            colvec <- c("grey10", "grey25", "grey40", "grey60", "grey80", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red")
                            matlab <- c(">0.9", 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005)
      
                            allcols[which(allcols1 >= 0.8)] <- "grey10"
                            allcols[which((allcols1  < 0.8) & (allcols1 >=0.7))] <- "grey25"
                            allcols[which((allcols1  < 0.7) & (allcols1 >=0.6))] <- "grey40"
                            allcols[which((allcols1  < 0.6) & (allcols1 >=0.5))] <- "grey60"
                            allcols[which((allcols1  < 0.5) & (allcols1 >=0.4))] <- "grey80"
                            allcols[which((allcols1  < 0.4) & (allcols1 >=0.3))] <- "SlateGray3"
                            allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "SlateGray1"
                            allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "DarkSeaGreen1"
                            allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "yellow"
                            allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "orange"
                            allcols[which((allcols1  < 0.01) & (allcols1 >=0.005))] <- "tomato"
                            allcols[which((allcols1  < 0.005))] <- "red"                                            
                      }
                      
                      if (col.scheme == "highsign") {
                      colvec <- c("grey45", "grey55", "SlateGray3", "SlateGray1", "DarkSeaGreen1", "yellow", "orange", "tomato", "red", "VioletRed1", "DeepPink", "DeepPink2")
                            matlab <- c(">0.3", 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001)
      
                            allcols[which(allcols1 >= 0.3)] <- "grey45"
                            allcols[which((allcols1  < 0.3) & (allcols1 >=0.2))] <- "grey55"
                            allcols[which((allcols1  < 0.2) & (allcols1 >=0.1))] <- "SlateGray3"
                            allcols[which((allcols1  < 0.1) & (allcols1 >=0.05))] <- "SlateGray1"
                            allcols[which((allcols1  < 0.05) & (allcols1 >=0.01))] <- "DarkSeaGreen1"
                            allcols[which((allcols1  < 0.01) & (allcols1 >=0.001))] <- "yellow"
                            allcols[which((allcols1  < 0.001) & (allcols1 >=0.0001))] <- "orange"
                            allcols[which((allcols1  < 0.0001) & (allcols1 >=0.00001))] <- "tomato"
                            allcols[which((allcols1  < 0.00001) & (allcols1 >=0.000001))] <- "red"
                            allcols[which((allcols1  < 0.000001) & (allcols1 >=0.0000001))] <- "VioletRed1"
                            allcols[which((allcols1  < 0.0000001) & (allcols1 >=0.00000001))] <- "DeepPink"
                            allcols[which((allcols1  < 0.000000001))] <- "DeepPink2"                                            
                      }                      
                      
                      
                      #pie(rep(1,length(allcols)), col=allcols)                      
                      
                      whcol <- which(as.character(sort(carps[[j]])) %in% as.character(carpP))
                      Pcols <- allcols[whcol]
                      #pie(rep(1,length(Pcols)), col=Pcols)
                      
                      ## um dopplete zu vermeiden ##
                      or <- order(carpP,decreasing=F)
                      r <-  vector(length=length(carpP))
                      r[or] <- c(1:length(carpP))
                      ##############################

                      #######################################
                      par(mar=c(2.5,0.1,2.5,30))
                      image(matrix(r,nrow=1), y=1:ncol(carpP), x=1:nrow(carpP), col=Pcols, axes=F, xlab="", ylab="", main=FName, cex.main=0.6)
                      axis(1, 1:nrow(carpP), labels="", las=2, cex.axis=1.5, lwd.ticks=0)
                      axis(4, 1:ncol(carpP), labels=rn, las=2, cex.axis=0.8, lwd.ticks=0)
                      abline(h=((1:ncol(carpP))+0.5), col="grey", lwd=0.1)
                      abline(v=((1:nrow(carpP))+0.5), col="grey", lwd=0.1)
                      box(lty = 'solid', col = 'grey3')
                      #######################################
      
                      ## draw z-Value colorscale on bottom ##
                      #par(mar=c(2, 10, 3.5, 35))
                      #matlab <- pretty(sort(as.vector(as.matrix(carp))),n=10)
                      #image(as.matrix(matlab), x=1:length(matlab),y=1, col=hcols, axes=F, xlab="", ylab="")
                      #box(lty = 'solid', col = 'black')
                      #axis(1, 1:length(matlab), labels=as.character(matlab), las=1, cex.axis=0.7, lwd.ticks=1)
                      #######################################
      
                      ## draw P-Value colorscale on bottom ##
                      par(mar=c(4, 0, 0, 25))  #par(mar=c(2, 2, 3.5, 3))
                      image(as.matrix(c(1:length(matlab))), x=1:length(matlab),y=1, col=colvec, axes=F, xlab="", ylab="")
                      box(lty = 'solid', col = 'black')
                      axis(1, 1:length(matlab), labels=as.character(matlab), las=3, cex.axis=0.8, lwd.ticks=1)
                       
                      #matlab <- pretty(sort(as.vector(as.matrix(carps[[j]]))),n=6)
                      #image(as.matrix(matlab), x=1:length(matlab),y=1, col=allcols, axes=F, xlab="", ylab="")
                      #box(lty = 'solid', col = 'black')
                      #axis(1, 1:length(matlab), labels=as.character(round(10^-matlab,3)), las=3, cex.axis=0.5, lwd.ticks=1)
                      #######################################
                      #graphics.off()
                 }
          }
      graphics.off()
     }
  if( Plot == "NO" ) { carpList }
}
###################################################



