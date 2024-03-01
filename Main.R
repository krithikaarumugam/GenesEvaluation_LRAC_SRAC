#KA on 23-07-2020
#!/usr/bin/env Rscript

#requires python3
#install.packages("argeparse")

library(argparse)
parser <- ArgumentParser()
#Parameters
#all mandatory

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="verbose [default]")
parser$add_argument("-q", "--quietly", action="store_false",dest="verbose", help="no verbose")
parser$add_argument("-Lfaa","--LRACfaa", type="character", help="Protein FASTA file from Prokka from LRAC", metavar="FASTA file")
parser$add_argument("-Sfaa","--SRACfaa", type="character", help="Protein FASTA file from Prokka from SRAC bin", metavar="FASTA file")
parser$add_argument("-d", "--diamondoutput", type="character", help="Diamond output file in tab format with 14 columns (lrac_ORFid, lrac_ORFlen, bin_srac_ORFid, bin_srac_ORFlen, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)", metavar="DIAMOND file")
parser$add_argument("-Lgff","--LRACgff", type="character", help="Protein FASTA file from Prokka (grep'ed for CDS)", metavar="gff file")
parser$add_argument("-Sgff","--SRACgff", type="character", help="Protein FASTA file from Prokka (grep'ed for CDS)", metavar="gff file")
parser$add_argument("-b","--blastn2.fil", type="character", help="blastn2 from srac2lrac filtered for al2ql>=0.95", metavar=".Rdata file")
args <- parser$parse_args()


if(is.null(args$LRACfaa) ==TRUE | is.null(args$SRACfaa) ==TRUE | is.null(args$diamondoutput) ==TRUE | is.null(args$LRACgff) ==TRUE | is.null(args$SRACgff) ==TRUE | is.null(args$blastn2.fil) == TRUE) 
{
	stop("Check arguments...!\n")
}


sink("genequality.log", append=TRUE, split=TRUE)  # for screen and log
cat("\n")
cat(paste(Sys.time()))
cat("

 **************************************************************************************************************************/      \\
                                                                                                                          |        |
                                                                                                                          |        |
                                                                                                                          |        |  
   GGGGGGGGG   EEEEEEEEE  NN      N  EEEEEEEEE     QQQQQQQQQ  U       U    AAAAA    L          I  TTTTTTTTT  Y       Y 	  \\        / 
   G           E          N N     N  E             Q       Q  U       U   A     A   L          I      T       Y     Y      \\      /  
   G           E          N  N    N  E             Q       Q  U       U  A       A  L          I      T         YYY         \\____/    
   G    GGGG   EEEEEEE    N   N   N  EEEEEEEE      Q    Q  Q  U       U  AAAAAAAAA  L          I      T         YYY         /    \\       
   G       G   E          N    N  N  E             Q     Q Q  U       U  A       A  L          I      T          Y         /      \\ 
   G       G   E          N     N N  E             Q      QQ  U       U  A       A  L          I      T          Y        /        \\
   GGGGGGGGG   EEEEEEEEE  N      NN  EEEEEEEEE     QQQQQQQQQ  UUUUUUUUU  A       A  LLLLLLLLL  I      T          Y        |        |        
                                                            Q                                                             |        |
                                                                                                                          |        |
 **************************************************************************************************************************\\      / \n")




#********************
cat(paste(Sys.time()," : Loading libraries... \nlibrary(ape)\nlibrary(dplyr)\nlibrary(scales)\n"))	
#********************

#following libraries are used
library(ape)
library(dplyr)
library(scales)
library(igraph)

#********************
cat(paste(Sys.time()," : Loaded libraries successfully!\n"))
cat(paste(Sys.time()," : You are here:",getwd(),"\n"))
#read diamond output and add column names. Query is LRACfaa and Subject is SRACfaa
cat(paste(Sys.time()," : Reading in diamond output of lracfaa and sracfaa\n"))
#********************

lracfaa_sracfaa_diamond<-read.table(args$diamondoutput,sep="\t",header=FALSE,stringsAsFactors=FALSE)
colnames(lracfaa_sracfaa_diamond)<-c("lrac_ORFid","lrac_ORFlen","bin_srac_ORFid","bin_srac_ORFlen","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
dim_lracfaa_sracfaa_diamond<-dim(lracfaa_sracfaa_diamond)

#********************
cat(paste(Sys.time()," : dim of lracfaa_sracfaa_diamond is:", dim_lracfaa_sracfaa_diamond[1], dim_lracfaa_sracfaa_diamond[2], "\n"))
cat(paste(Sys.time()," : Length of lracfaa_sracfaa_diamond$lrac_ORFid:", length(lracfaa_sracfaa_diamond$lrac_ORFid),"\n"))
cat(paste(Sys.time()," : Number of LRAC genes having a match to a gene in SRAC:",length(unique(lracfaa_sracfaa_diamond$lrac_ORFid)),"\n"))
#********************

if(length(lracfaa_sracfaa_diamond$lrac_ORFid) == length(unique(lracfaa_sracfaa_diamond$lrac_ORFid)))
  {
	   cat(paste(Sys.time()," : No duplicates in lrac_ORFid!\n"))
  }else
  {
	   cat(paste(Sys.time()," : There's a duplicate element in lrac_ORFid!! Check lracfaa_sracfaa_diamond$lrac_ORFid\n"))
  }

#********************
cat(paste(Sys.time()," : Reading LRAC.faa and SRAC.faa files from Prokka...\n"))
#********************

lracfaa<-read.FASTA(args$LRACfaa,type="AA")
sracfaa<-read.FASTA(args$SRACfaa,type="AA")

#********************
cat(paste(Sys.time()," : Number of LRAC genes predicted by prokka:", length(lracfaa),"\n"))
cat(paste(Sys.time()," : Number of SRAC genes predicted by prokka:", length(sracfaa), "\n"))
cat(paste(Sys.time()," : Number of SRAC genes matched to an LRAC gene:", length(unique(lracfaa_sracfaa_diamond$bin_srac_ORFid)),"\n"))
cat(paste(Sys.time()," : Number of SRAC genes not mapped to LRACgenes:", length(sracfaa) - length(unique(lracfaa_sracfaa_diamond$bin_srac_ORFid)),"\n"))
cat(paste(Sys.time()," : Number of LRAC genes not mapped to SRACgenes:", length(lracfaa) - length(unique(lracfaa_sracfaa_diamond$lrac_ORFid)),"\n"))
cat(paste(Sys.time()," : Reading LRAC gff file from Prokka...\n"))
#********************

#read in LRAC gff files after filtering for CDS lines
LRAC_gff<-read.gff(args$LRACgff)

#********************
cat(paste(Sys.time()," : Reading SRAC gff file from Prokka...\n"))
#********************

#read in SRAC gff files after filtering for CDS lines
SRACbin_gff<-read.gff(args$SRACgff)

#********************
cat(paste(Sys.time()," : Loading blastn2.fil.RData from srac2lrac...\n"))
#********************

#load blastn2.fil from srac2lrac
load(args$blastn2.fil)

#********************
cat(paste(Sys.time()," : Read all 6 input files successfully! CALLING GENE QUALITY FUNCTION now...\n"))
#********************


source(file="Genequality.R")

GeneQuality_LRAC_SRAC(lracfaa,sracfaa,lracfaa_sracfaa_diamond,LRAC_gff,SRACbin_gff,blastn2.uni.fil)

sink()


