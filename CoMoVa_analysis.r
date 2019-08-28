################################################################################
#
#
#


library(gplots)
library(seqLogo)

print("Loading orthologs...")
if(!exists("orthologs")) {
	load("./orthologDB/orthologs.Rworkspace")
}

print("Loading CoMoVa...")

source("CoMoVa.r")

print("Loading Annotations...")
if(!exists("annotation_s")) {
  annotation_s = as.matrix(read.csv("tair_annotation.csv", as.is=TRUE))[,2:3]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,2)!= "At", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "NAC", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "PUB", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "ECP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATCOR", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "ATHB", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "AXR5", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "ZIP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,7)!= "ATAUX2-", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATIRT", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATCTL", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,6)!= "ATEXPA", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,6)!= "ATHEXP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "EXP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,6)!= "ATNCED", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATCAO", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "SIS", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "CH1", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "STO1", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "ATRR", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "ASL", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "IBC", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "MEE", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "TTR", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,2)!= "RR", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATTLL", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATLUP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "ATHH", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATFRO", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "GAMMA", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "SAMDC", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "HMGB", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "G-H2A", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "WSIP", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "H2AX", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATERF", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "ATAF1", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,5)!= "TMAC2", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "AXR3", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,6)!= "GATA18", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "APRR", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,4)!= "AIP1", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "HON", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "SHD", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "TL1", ]
  annotation_s = annotation_s[ substr(annotation_s[,2],1,3)!= "MNP", ]
  
  annotation_s[,2] = toupper(annotation_s[,2])

  annotation_s = sapply(unique(annotation_s[,"TAIR"]), function(x) { paste(unique(head(sort(annotation_s[annotation_s[,"TAIR"]==x,"SYMBOL"]),n=2)), collapse="\\",sep="")})
  names(annotation_s) = sapply(names(annotation_s), function(x) { ifelse(regexpr("\\.",x)[[1]]>0, substr(x,1,regexpr("\\.",x)-1),x) })
  annotation_s["AT3G55580"]="TCF1"
  annotation_s["AT4G16660"]="HSP70"
  annotation_s["AT1G49780"]="PUB26"
  annotation_s["AT3G53830"]="RCC1"
  annotation_s["AT5G44550"]="CASPL1B1"
  annotation_s["AT1G07820"]="HIS4Like"
  annotation_s["AT1G07660"]="HIS4Like"
  annotation_s["AT3G45930"]="HIS4Like"
  annotation_s["AT1G06760"]="HIS1.1"
  annotation_s["AT4G37180"]="UIF1"
  annotation_s["AT2G30620"]="HIS1.2"
  annotation_s["AT4G37390"]="GH3.2"
  annotation_s["AT2G23830"]="PVA31"
  annotation_s["AT1G60190"]="ATPUB19"
  annotation_s["AT1G19210"]="ERF017"
  annotation_s["AT1G69570"]="CDF5" 
}


PlotConHM <- function(pval, cutoff=30, minvalue = 10, maxval=30, trim= FALSE, returnMatrix=FALSE, cexCol = 1.8, cexRow = 1, annotate=FALSE, plotHM=TRUE) {

  by=colorRampPalette(rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))) 

  selected = rownames(pval)[apply(pval,1,max)>= cutoff]
  pval_sel = pval[selected,]

  colsep=c(0,ncol(pval))
  if(ncol(pval)<30) {
    colsep = 0:ncol(pval)
  }
  rowsep=c(0,length(selected))
  if(length(selected)<50) {
    rowsep=0:length(selected)
  }
  
  if(trim) {
    pval_sel = pval_sel[, apply(pval_sel, 2, max)>=cutoff]
    order_col = intersect(order_col, colnames(pval_sel))
  }  
  pvalfil = pval_sel
  pvalfil[pvalfil<cutoff]=0
  corder = order(apply(pval_sel,2,max), decreasing = TRUE)
  if(annotate) { 
     rownames(pval_sel) = sapply(rownames(pval_sel), function(x) { ifelse(x %in% names(annotation_s), annotation_s[x], x)})
  }
  roworder = rev(do.call(order, as.data.frame(pvalfil[,corder])))

  if(plotHM) {
     heatmap.2(pval_sel[roworder,corder], col=c("#FFFFFF",by(19)), scale="none", Rowv=NULL, Colv=NULL, density.info = "none", trace="none", breaks=seq(minvalue,min(max(pval),maxval),(min(max(pval),maxval)-minvalue)/20),
               margins=c(10,10), dendrogram = "none", rowsep = rowsep, colsep = colsep, sepcolor = "black",
               sepwidth = c(0.0001,0.0001), keysize = 1, key.title = "", key.xlab = "-log(Pval)", cexCol= cexCol, cexRow = cexRow, adjCol= c(0.9,0.4))
   }
  if(returnMatrix) {
     pval_sel
  }
}

getSuggestedCutoff <- function(background, backgroundcutoff=1) {
  fd=vector()
  for(i in 1:100) { fd[i] = sum(background >=i ) }
  return(min(which(fd<= backgroundcutoff)))
}

PlotPWMForMotif <- function(redb, locus, variant , orthologlist  = c(monocots, dicots)) {
  orthlist = getAllOrthologsByRE(locus, orthologs, redb, orthologlist=orthologlist)
  tab = do.call(rbind,lapply(names(orthlist), function(x) { redb[[x]][redb[[x]][,1]== orthlist[x] & redb[[x]][,"Variant"]==variant,]}))
  
  if(nrow(tab) <= 1) {
    return(0)
  }
  #if there are two motifs, prefer the last.
  
  # reverse if needed
  tab[tab[,"Rev"]=="R", "Region"] = sapply(tab[tab[,"Rev"]=="R","Region"], function(x) 
    { as.character(reverseComplement(DNAString(x))) })
  # drop no info
  tab = tab[regexpr("-",tab[,"Region"])== -1,]
  # start by picking the first one, if there are multiple ones. We will optimize later
  tab = tab[order(tab[,"Pos"], decreasing = FALSE),]

  #select promoters with multiple repeats
  multVarLoci = unique(tab$Locus[duplicated(tab$Locus)])
  
  # start with the default set - just the first one
  bestset = !duplicated(tab$Locus)
  curset=bestset
  curcm = consensusMatrix(DNAStringSet(tab[curset,"Region"]), as.prob = TRUE)[c("A","C","G","T"),]
  bestcm = curcm
  #try multiple ones
  for(curLocus in multVarLoci) {
    possible = which(tab$Locus== curLocus)
    first = possible[1]
    for(curpos in possible[2:length(possible)]) {
      #try to replace.
      curset[possible]=FALSE
      curset[curpos]=TRUE
      curcm = consensusMatrix(DNAStringSet(tab[curset,"Region"]), as.prob = TRUE)[c("A","C","G","T"),]
      if(sum(apply(curcm,2,max))>sum(apply(bestcm,2,max))) {
       bestset = curset
        bestcm = curcm
      } else {
        curset=bestset
      }
    }
  }
  positions = tab[bestset,"Pos"]
  cm = consensusMatrix(DNAStringSet(tab[bestset,"Region"]), as.prob = TRUE)[c("A","C","G","T"),]

  #normalize  
  for(i in 1:ncol(cm)) { cm[,i] = cm[,i]/sum(cm[,i])}

  #makePwm
  pwm_gene=makePWM ( cm)
  print(seqLogo(pwm_gene, xaxis = FALSE, yaxis=FALSE))
  
  return(list(pwm=pwm_gene, num=sum(bestset)))
}


