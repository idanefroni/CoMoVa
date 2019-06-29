
if(file.exists("orthologDB/orthologs.Rworkspace")) {
        cat("Loading ortholog database...\n")
        load("orthologDB/orthologs.Rworkspace")
} else {
	stop("ERROR: Can't find ortholog file")
}
suppressMessages(source("CoMoVa.r"))

genelist = sort(rownames(orthologs[["Athaliana"]]))
args = commandArgs(trailingOnly=TRUE)

if(length(args)<4) {
	stop("ERROR: missing paramters. Syntax: buildConservationMatrix <motifname> <promoterlength> <startgene> <endgene>")
} else {
	motifname = args[1]
	promoterlength = as.numeric(args[2])
	a = as.numeric(args[3])
	b = as.numeric(args[4])
}

dumpConservation <- function(redb, genes, outFile, header) {
  dotcounter=0
  
  brafile = paste("tmp/bra_", outFile,sep="")
  dicotsfile = paste("tmp/dicots_", outFile,sep="")
  monocotsfile= paste("tmp/monocots_",outFile,sep="")
  angiofile= paste("tmp/angio_",outFile,sep="")

  if(header) {
    cat(c("Locus", names(getConservationData(genes[1], redb)$conserv)), sep=",", file=angiofile, append=TRUE)
    cat("\n", file=angiofile, append=TRUE)
    cat(c("Locus", names(getConservationData(genes[1], redb)$conserv)), sep=",", file=dicotsfile, append=TRUE)
    cat("\n", file=dicotsfile, append=TRUE)
    cat(c("Locus", names(getConservationData(genes[1], redb)$conserv)), sep=",", file=monocotsfile, append=TRUE)
    cat("\n", file=monocotsfile, append=TRUE)
    cat(c("Locus", names(getConservationData(genes[1], redb)$conserv)), sep=",", file=brafile, append=TRUE)
    cat("\n", file=brafile, append=TRUE)   
  }
  tmp=sapply(genes, function(x) {
    cv = getConservationData(x,redb)
    cat(".")
    dotcounter=dotcounter+1
    if(dotcounter==100) {
    	cat("\n")
    	dotcounter=0
    }
    
    cat(c(x,cv$conserv),sep=",",  file=angiofile, append=TRUE)
    cat(c(x,cv$dicots_conserv),sep=",",  file=dicotsfile, append=TRUE)
    cat(c(x,cv$monocots_conserv),sep=",",  file=monocotsfile, append=TRUE)
    cat(c(x,cv$brassica_conserv),sep=",",  file=brafile, append=TRUE)
    
    cat("\n",  file=angiofile, append=TRUE)
    cat("\n",  file=dicotsfile, append=TRUE)
    cat("\n",  file=monocotsfile, append=TRUE)
    cat("\n",  file=brafile, append=TRUE)
    })
}


# load data

smalldatafile = paste0("REDB/", motifname , ".small")

if(file.exists(smalldatafile)) {
	datafile=smalldatafile
} else {
	datafile = paste0("REDB/", motifname , ".data")
}
	
if(file.exists(datafile)) {
    cat("Loading RE database...")
    load(file=datafile)
    redb=get(motifname)
    rm(list=motifname)
    cat("Done.\n")
} else {
	stop("ERROR: Can't find RE database file")
}

# filter by promoter length and drop memory-greedy fields
cat(paste0("Filter by promoter length [" , promoterlength, "]...\n"))
target_genelist = genelist[a:(min(c(length(genelist),b)))]

for(i in names(redb)) {
	redb[[i]]= redb[[i]][redb[[i]][,"Pos"]< promoterlength, c("Locus","Variant","Pos","Rev")]
}
gc(verbose=FALSE, full=TRUE)

cat(paste0("Generating motif conservation data for ", length(target_genelist), " genes [" , a , "-", b,"]..."))
dumpConservation(redb, target_genelist, paste0(motifname, "_", promoterlength , "_", a, "_",b,".csv"), (a==1))
