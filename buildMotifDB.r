suppressMessages(source("CoMoVa.r"))

if(file.exists("orthologDB/orthologs.Rworkspace")) {
        cat("Loading ortholog database...\n")
        load("orthologDB/orthologs.Rworkspace")
} else {
	stop("ERROR: Can't find ortholog file")
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
	stop("ERROR: Must provide motif name")
} else {
	motifname = trimws(args[1])
}

loadREDB <- function(name, sourceForBackground=NULL, trim=TRUE) {
  reDB <- list()
  flist <- dir(path=paste0("./REDB/", name), pattern = "*.txt")
  cat("Processing...")
  for(af in flist) {
    sp = substr(af,1,which(strsplit(af,'')[[1]]=="_")[[1]]-1)
    cat(paste0(sp,","))
    reDB[[sp]] = read.csv(paste0("./REDB/",name,"/", af), header=TRUE, as.is=TRUE, row.names=NULL)
    if(trim) {
      reDB[[sp]] = reDB[[sp]][,c("Locus","Variant","Pos","Rev")]
    }
    rownames(reDB[[sp]]) = c()

#    if(sp=="Slycopersicum") {
#      reDB[[sp]][,1] = stripDot(reDB[[sp]][,1])
#    }
    
    reDB[[sp]] = unique(reDB[[sp]])
    
    #if this is a background database, downsample to original size
    if(!is.null(sourceForBackground)) {
      reDB[[sp]] = reDB[[sp]][sample(1:nrow(reDB[[sp]]), nrow(sourceForBackground[[sp]]), replace = TRUE  ),]
    }
  }
  cat("...Done.\n")
  reDB
}

reDB = loadREDB(motifname, trim=FALSE)
assign(motifname, reDB)
assign(paste0(motifname, "background") , loadREDB(paste0(motifname, "background"),sourceForBackground= reDB))

save(list=motifname, file=paste0("REDB/", motifname, ".data"))
save(list= paste0(motifname, "background"), file=paste0("./REDB/", motifname, "background.data"))
	
#make small version
reDB = loadREDB(motifname, trim=TRUE)
assign(motifname, reDB)
save(list=motifname, file=paste0("REDB/", motifname, ".small"))
