# CoMoVa

CoMoVa is a utility for identification of conserved DNA motif variants in gene promoters

################################################################################################
# Setting up

Obtain CoMoVa from https://github.com/idanefroni/CoMoVa

run setup using:

chmod +x setupCoMoVa
./setupCoMoVa

This will build the directory structure and setup the required permissions

To use CoMoVa you will need to setup the genome database and the response element database.

#####################################
#### Setting up the genome database

For each genome CoMoVa requires three files: a genome fasta file (gzipped), a GFF3 genome annotation file and a protein sequence file.
The three files must be located in the genome directory and follow the following naming convention (same as phytozome):
Genome file: <GenomeName>.fa.gz
GFF3 file: <Genomename>.gene.gff3
Protein file: <GenomeName>.protein_primaryTranscriptOnly.fa

For example, for the Genome Athaliana_167_TAIR10 the files will be:
Athaliana_167_TAIR10.fa.gz
Athaliana_167_TAIR10.gene.gff3
Athaliana_167_TAIR10.protein_primaryTranscriptOnly.fa

The github distribution already has the Arabidopsis genome in the directory.

Following obtaining the files, the genome has to be processed and orthology determined. This is done with the processGenome script

./processGenome <GenomeName>

Depending on your computing power, this can take a couple of minutes to an hour to complete for each genome. You can also run it as 
a batch for multiple genome using the processAllGenomes script.

The script will generate for each genome, a blast database in the blastdb directory, a promoter sequence file in the genome directory names
<Genomename>.upstream.fa, and a reciprocal blast tables in the orthologDB directory. It will also add generate an orthologs.Rworkspace in
the orthologDB directory that will be used by downstream analysis.

###########################################
##### Setting up the motif variant database

Motifs are located in the motifs directory and are defined in FASTA format. Each motif should have two FASTA files for the motif code and
for its background (see the paper for details).
The motif uses the IUPAC code for nucleotides with the following modification: N stands for nucleotide to be skipped. 
X stand for any base that should form the motif variant. See the example files and the paper for details.
The names for the files should be:
<Motifname>.fasta
<Motifname>background.fasta

Motif database is generated using the processMotif script:
./processMotif <Motifname>


###############################################
#### Computing motif conservation

Once the genomes and motif database is setup, we calculate motif conservation:

./CoMoVa <Motifname> <promoter length> >logs/<Motifname>.log

This script may take a considerable amount of time which depends on how common is the motif. Computing conservation for the auxRE
for 45 genomes using a 16 core Xeon 2.4GHz took ~24 hours. This will generate a conservation database for the response element. 


##############################################################
##### Data analysis

The scripts for data analysis are written in R. All basic functions are in the CoMoVa.r file.

Sample analysis R code:

# Loads CoMoVa and basic databases
source("CoMoVa_analysis.r")


# loading a particular RE database (replace "auxRE" with the particular RE in interest)
# first, load the background database and the motif conservation database

auxbackground <- loadcontable("auxREbackground", length=1000)
aux <- loadcontable("auxRE", length=1000, background=auxREbackground)

# To get the significantly conserved genes, first, determine the appropriate cutoff
#  the cutoff is determined based on the specific background distribution. It can (and should) be twicked manually.

cutoff = getSuggestedCutoff(auxbackground$angio_pval)

# Then, plot the conservation heatmap (similar to Fig 2A in the paper)

aux_cons_matrix = PlotConHM(aux$angio_pval, cutoff=cutoff, returnMatrix=TRUE, plotHM=TRUE)

# aux_cons_matrix will contain a matrix, where the rows are genes and columns motif variants. Values are the conservation p-values.

# Alternative analysis options
#
# * To analyze the raw conservation scores (instead of pval) use aux$angio instead of aux$angio_pval
#
# * To analysis just a subset of the tree use aux$dicots_pval, aux$monocots_pval or bra_pval (for Brasicacea). 
#
# * To return an annotated table having gene names (instead of locus ID), use annotate=TRUE in PlotConHM
#



# To plot the PWM for particular motifs (Fig 5):
# First, load the raw variant data

load("REDB\auxRE.data")

# Then compute PWM for a give locus (example is for IAA2)

PlotPWMforMotif(redb= auxRE, locus = "AT3G23030", motif="CC")





###############################################
### Version tracking

## V0.01 - First upload
## V0.02 - Added analysis functions




