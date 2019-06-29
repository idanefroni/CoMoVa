suppressMessages(source("CoMoVa.r"))

if(file.exists("orthologDB/orthologs.Rworkspace")) {
	cat("Loading ortholog database...\n")
	load("orthologDB/orthologs.Rworkspace")
}
keepBlastTable = FALSE
touched=FALSE

if(!exists("orthologs")) {
	orthologs=list()
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
	genome = ""
} else {
	genome = args[1]
}

bn <- dir(path="orthologDB", pattern= paste0(".*", genome, ".*txt$"))
species = substr(bn[grep("^ara2", bn)],5,100)
a = paste("ara2", species, sep="")
b= paste(substr(species, 1, nchar(species)-4), "2ara.txt",sep="")
speciesShort = substr(species,1,which(strsplit(species,'')[[1]]=="_")[[1]]-1)
names(speciesShort) = species

cat(paste0("Processing orthologs for ", speciesShort , "\n"))

if(!(speciesShort %in% names(orthologs))) {
	if(!exists("blasttables")) {
		blasttables=list()
	}

	for(i in bn) {
		cat(paste0("Load table ", i ,"..."))
		if(!(i %in% names(blasttables))) {
			blasttables [[i]] = eval2ranks(read.table(paste("orthologDB//", i, sep=""), header=FALSE, as.is=TRUE))
			cat("Done\n")
		} else {
			cat("already loaded, skipped.\n") 
		}
	}
	cat(paste0("Building orthologs for ", speciesShort , "..."))
	orthologs[[speciesShort]] = buildOrthologs(as.matrix(blasttables[[a]]),as.matrix(blasttables[[b]]))
	cat("Done\n")
	touched=TRUE

} else {
	cat(paste0("Ortholog table " , speciesShort, " is already loaded, skipped...\n")) 
}

suppressWarnings( rm(bn,a,b,i))

if(!keepBlastTable) {
	suppressWarnings( rm(blasttables))
}
if(touched) {
	cat("Saving...")
	save.image("orthologDB/orthologs.Rworkspace")
	cat("Done\n")
}
