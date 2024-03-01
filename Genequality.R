#you need to know LRAC id in blastn2.fil - user input
#Gene Quality Function
GeneQuality_LRAC_SRAC<-function(Lracfaa,Sracfaa,LracfaaSracfaa.diamond,LRACgff,SRACgff,blastn2.fil)
{
#********************
cat(paste(Sys.time()," : Editing attributes column in gff files and adding extra columns to include ORFid and annotation...\n"))
#********************	
#LRAC
#extracting ORF id from the attributes column and create a new column named 'ID'
LRACgff<-cbind(LRACgff, ID=sapply(strsplit(LRACgff$attributes,";"),FUN=function(x) {x[1]}),stringsAsFactors = FALSE)
#splitting 'ID' to extract the ORF id and create a new column named 'ID.cut'
LRACgff<-cbind(LRACgff, ID.cut=sapply(strsplit(LRACgff$ID,"="),FUN=function(x){x[2]}),stringsAsFactors = FALSE)
#include product
LRACgff<-cbind(LRACgff, product=sapply(strsplit(LRACgff$attributes,"product="),FUN=function(x){x[2]}),stringsAsFactors = FALSE)
#SRAC
#extracting ORF id from the attributes column and create a new column named 'ID'
SRACgff<-cbind(SRACgff, ID=sapply(strsplit(SRACgff$attributes,";"),FUN=function(x) {x[1]}),stringsAsFactors = FALSE)
#splitting 'ID' to extract the ORF id and create a new column named 'ID.cut'
SRACgff<-cbind(SRACgff, ID.cut=sapply(strsplit(SRACgff$ID,"="),FUN=function(x){x[2]}),stringsAsFactors = FALSE)
#include product
SRACgff<-cbind(SRACgff, product=sapply(strsplit(SRACgff$attributes,"product="),FUN=function(x){x[2]}),stringsAsFactors = FALSE)

#********************
cat(paste(Sys.time()," : Merging diamond output and GFF files of SRAC and LRAC...\n"))
#********************

#merge LRACgff with lracfaa_sracfaa_diamond
merge_lracfaa_sracfaa_diamond_LRACgff<-merge(LracfaaSracfaa.diamond,LRACgff,by.x="lrac_ORFid",by.y="ID.cut")#all.y=TRUE)

#********************
cat(paste(Sys.time()," : dim of merge_lracfaa_sracfaa_diamond_LRACgff is:", dim(merge_lracfaa_sracfaa_diamond_LRACgff)[1], "\n") )
#********************

#merge above and SRACgff based on sracorfid 
merge_lracfaa_sracfaa_diamond_LRACgff_SRACgff<-merge(merge_lracfaa_sracfaa_diamond_LRACgff,SRACgff,by.x="bin_srac_ORFid",by.y="ID.cut")

#********************
cat(paste(Sys.time()," : dim of merge_lracfaa_sracfaa_diamond_LRACgff_SRACgff is:", dim(merge_lracfaa_sracfaa_diamond_LRACgff_SRACgff)[1], "\n") )
cat(paste(Sys.time()," : check if dim(merge_lracfaa_sracfaa_diamond_LRACgff_SRACgff)[1] == dim(merge_lracfaa_sracfaa_diamond_LRACgff)[1]\n"))
cat(paste(Sys.time()," : delete columns not used to increase readability of df","\n"))
#********************

df = subset(merge_lracfaa_sracfaa_diamond_LRACgff_SRACgff,select = -c(source.x,type.x,score.x,phase.x,attributes.x,source.y,type.y,score.y,phase.y,attributes.y),stringsAsFactors = FALSE)
#change class, factor -> character
df$seqid.y<-as.character(df$seqid.y)

#********************
cat(paste(Sys.time()," : everything '.x' is long read, everything '.y' refers to an SRAC", "\n"))
cat(paste(Sys.time()," : head(df)\n"))
print(head(df))
cat(paste(Sys.time()," : subset based on \"1\", since \"1\" is chr\n"))
#********************

blastn2.uni.fil.1<-subset(blastn2.fil, sseqid=='1')

#********************
cat(paste(Sys.time()," : modify df$seqid.y to match qseqid in blastn2\n"))
#********************

df<-cbind(df,shortseq=sapply(strsplit(df$seqid.y,"_"),FUN=function(x){paste(x[[1]],x[[2]],sep="_")}))

#********************
cat(paste(Sys.time()," : head(df)\n"))
print(head(df))
cat(paste(Sys.time()," : merge df and blastn2.uni.fil.1. all.x is set to true to maintain all rows of the original diamond output\n"))
#********************

merge_df_blastn2.uni.fil.1<-merge(df,blastn2.uni.fil.1,by.x="shortseq",by.y="qseqid",all.x=TRUE)

#********************
cat(paste(Sys.time()," : check for cases where SRAC do not align to LR-chr using blastn2.uni.fil results: table(is.na(merge_df_blastn2.uni.fil.1$sstart.y))","\n"))
print(table(is.na(merge_df_blastn2.uni.fil.1$sstart.y)))
cat(paste(Sys.time()," : SRAC sequences not aligning to LR-chr: data.frame(table(merge_df_blastn2.uni.fil.1[is.na(merge_df_blastn2.uni.fil.1$sstart.y),\"seqid.y\"]))","\n"))
print(data.frame(table(merge_df_blastn2.uni.fil.1[is.na(merge_df_blastn2.uni.fil.1$sstart.y),"seqid.y"])))
cat(paste(Sys.time()," : Calculating midpoints..\n"))
cat(paste(Sys.time()," : calculate midpoints of genes on the LR-chr\n"))
#********************

merge_df_blastn2.uni.fil.1$mplracfaa<-(merge_df_blastn2.uni.fil.1$start.x + merge_df_blastn2.uni.fil.1$end.x)/2

#********************
cat(paste(Sys.time()," : calculate mid-points of SRAC-genes mapping to LR-chr\n"))
#********************

merge_df_blastn2.uni.fil.1$mpsracfaa.all<-((merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$end.y))/2

#********************
cat(paste(Sys.time()," : calculate mid-points of SRAC-genes mapping to LR-chr considering strand information but not order of genes \n"))
#********************

merge_df_blastn2.uni.fil.1$mpsracfaa.strand<-with(merge_df_blastn2.uni.fil.1,ifelse(merge_df_blastn2.uni.fil.1$sstart.y>merge_df_blastn2.uni.fil.1$send.y,(((merge_df_blastn2.uni.fil.1$send.y+merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$send.y+merge_df_blastn2.uni.fil.1$end.y))/2),(((merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$end.y))/2)))
#--this is the approach we tried on Tuesday afternoon (30 June)
#--for SRAC aligned to -ve strand of LR-chr, get the correct positions by subtracting the start/end of a gene from strt position of the (aligned) contig on the LR-chr
#mptest<-with(merge_df_blastn2.uni.fil.1,ifelse(merge_df_blastn2.uni.fil.1$sstart.y>merge_df_blastn2.uni.fil.1$send.y,(((merge_df_blastn2.uni.fil.1$sstart.y-merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y-merge_df_blastn2.uni.fil.1$end.y))/2),(((merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$end.y))/2)))
#********************
cat(paste(Sys.time()," : calculate mid-points of SRAC-genes mapping to LR-chr considering strand information and order of genes on the strand\n"))
#********************

merge_df_blastn2.uni.fil.1$mpsracfaa.strandnew<-with(merge_df_blastn2.uni.fil.1,ifelse(merge_df_blastn2.uni.fil.1$sstart.y>merge_df_blastn2.uni.fil.1$send.y,(((merge_df_blastn2.uni.fil.1$sstart.y-merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y-merge_df_blastn2.uni.fil.1$end.y))/2),(((merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$start.y)+(merge_df_blastn2.uni.fil.1$sstart.y+merge_df_blastn2.uni.fil.1$end.y))/2)))


#********************
cat(paste(Sys.time()," : CHECK IF GENES ARE RUNNING OFF ENDS OF CONTIGS - start of SRACfaa must be > 1 and end of SRACfaa must be less than length of SRAC. Getting length of SRAC from seqid.y since qlen is true only if present in blastn2.uni :-o, so create a new column seqid.y.length\n"))
#********************

merge_df_blastn2.uni.fil.1<-cbind(merge_df_blastn2.uni.fil.1, seqid.y.length=sapply(strsplit(merge_df_blastn2.uni.fil.1$seqid.y,"_"),FUN=function(x) {x[4]}),stringsAsFactors = FALSE)
merge_df_blastn2.uni.fil.1$seqid.y.length<-as.numeric(merge_df_blastn2.uni.fil.1$seqid.y.length)
merge_df_blastn2.uni.fil.1$generunoffSRAC<- ifelse(merge_df_blastn2.uni.fil.1$start.y > 1 & merge_df_blastn2.uni.fil.1$end.y < merge_df_blastn2.uni.fil.1$seqid.y.length,'TRUE','FALSE')

#********************
cat(paste(Sys.time()," : Number of genes running off the edge: ", table(merge_df_blastn2.uni.fil.1$generunoffSRAC)[1], ",\t", "Number of genes not running off the edge: ", table(merge_df_blastn2.uni.fil.1$generunoffSRAC)[2],"\n"))
cat(paste(Sys.time()," : Index of genes running off the edge of contig: ", which(grepl("FALSE",merge_df_blastn2.uni.fil.1$generunoffSRAC)),"\n"))
cat(paste(Sys.time()," : Number of cases where LRAC gene annotation doesn't match its SRAC gene annotation: ", table(merge_df_blastn2.uni.fil.1$product.x == merge_df_blastn2.uni.fil.1$product.y)[1],",\t", "Number of cases where LRAC gene annotation matches its SRAC gene annotation", table(merge_df_blastn2.uni.fil.1$product.x == merge_df_blastn2.uni.fil.1$product.y)[2],"\n"))
cat(paste(Sys.time()," : CHECK GENES WHERE QLEN/SLEN == 1\n"))
#********************


merge_df_blastn2.uni.fil.1.watsonratio.eq1<-subset(merge_df_blastn2.uni.fil.1,(merge_df_blastn2.uni.fil.1$lrac_ORFlen/merge_df_blastn2.uni.fil.1$bin_srac_ORFlen) == 1)

#********************
cat(paste(Sys.time()," : Number of LRAC genes having a match to SRAC genes with genelengthratio equal to 1: ", dim(merge_df_blastn2.uni.fil.1.watsonratio.eq1)[1],"\n"))
cat(paste(Sys.time()," : Number of LRAC genes having a match to SRAC genes with genelengthratio not equal to 1: ", dim(merge_df_blastn2.uni.fil.1)[1] -  dim(merge_df_blastn2.uni.fil.1.watsonratio.eq1)[1],"\n"))
cat(paste(Sys.time()," : CHECK GENES WHERE QLEN/SLEN > 1\n"))
#********************

merge_df_blastn2.uni.fil.1.watsonratio.gt1<-subset(merge_df_blastn2.uni.fil.1,(merge_df_blastn2.uni.fil.1$lrac_ORFlen/merge_df_blastn2.uni.fil.1$bin_srac_ORFlen) > 1)

#********************
cat(paste(Sys.time()," : Number of LRAC genes having a match to SRAC genes with genelengthratio > 1: ", dim(merge_df_blastn2.uni.fil.1.watsonratio.gt1)[1],"\n"))
cat(paste(Sys.time()," : CHECK GENES WHERE QLEN/SLEN > 2\n"))
#********************


merge_df_blastn2.uni.fil.1.watsonratio.gteq2<-subset(merge_df_blastn2.uni.fil.1,(merge_df_blastn2.uni.fil.1$lrac_ORFlen/merge_df_blastn2.uni.fil.1$bin_srac_ORFlen) >= 2)

#********************
cat(paste(Sys.time()," : Number of LRAC genes having a match to SRAC genes with genelengthratio > 2: ", dim(merge_df_blastn2.uni.fil.1.watsonratio.gteq2)[1],"\n"))
#cat(paste(Sys.time()," : CHECKING FOR Glg genes..\n"))
#********************

#merge_df_blastn2.uni.fil.1.glginproduct.x<-subset(merge_df_blastn2.uni.fil.1, grepl("Glg", merge_df_blastn2.uni.fil.1$product.x))

#********************
#cat(paste(Sys.time()," : Number of Glg genes: ",dim(merge_df_blastn2.uni.fil.1.glginproduct.x)[1],"\n"))
cat(paste(Sys.time()," : TESTING FOR SPLIT GENES; start with identifying duplicates in bin_srac_ORFid, meaning multiple lrac_ORFid could map to same bin_srac_ORFid? (but they don't necessarily mean split genes so check sracfaa co-ordinates)\n"))
#********************

merge_df_blastn2.uni.fil.1$bin_srac_ORFid<-as.character(merge_df_blastn2.uni.fil.1$bin_srac_ORFid)
#using dplyr to include #bin1_srac_ORFid >1
merge_df_blastn2.uni.fil.1.subset_dup.sracfaa <- merge_df_blastn2.uni.fil.1 %>% group_by(bin_srac_ORFid) %>% filter(n()>1) %>%as.data.frame() 

#********************
cat(paste(Sys.time()," : dim[1] of dataframe subset (or subsetted?is a word?) for multiple occurence of srac)",dim(merge_df_blastn2.uni.fil.1.subset_dup.sracfaa)[1],"\n"))
cat(paste(Sys.time()," : subsetting merge_df_blastn2.uni.fil.1.subset_dup.sracfaa where mpsracfaa.strandnew != NA ...\n"))
#********************

tmp<-merge_df_blastn2.uni.fil.1.subset_dup.sracfaa %>% filter(!is.na(mpsracfaa.strandnew)) %>% as.data.frame()

#********************
cat(paste(Sys.time()," : Number of SRAC genes occurring multiple times: ",length(unique(tmp$bin_srac_ORFid)),"\n"))
cat(paste(Sys.time()," : Subsetting LRAC genes without a match to SRAC genes in diamond output from lracgff file...\n"))
#********************

no_diamondmatch<-subset(LRACgff, !(LRACgff$ID.cut %in% LracfaaSracfaa.diamond$lrac_ORFid))

#********************
cat(paste(Sys.time()," : Number of LRAC genes not mapping to SRAC genes: ", dim(no_diamondmatch)[1], "\n"))
cat(paste(Sys.time()," : QUANTILES\n")) 
cat(paste(Sys.time()," : Adding the difference in midpoint positions as a column...\n"))
#********************

merge_df_blastn2.uni.fil.1$diff<-abs(merge_df_blastn2.uni.fil.1$mplracfaa-merge_df_blastn2.uni.fil.1$mpsracfaa.strandnew)

#********************
cat(paste(Sys.time()," : Adding protein length ratio (plr) data as another column...\n"))
#********************

merge_df_blastn2.uni.fil.1$plr<-merge_df_blastn2.uni.fil.1$lrac_ORFlen/merge_df_blastn2.uni.fil.1$bin_srac_ORFlen

#********************
cat(paste(Sys.time()," : if we use absolute value of the difference in positions, we only need to tag aligned proteins in the top end of the distribution, not bottom end, so we only need two colours. Take a look at every half percentile difference to get a more fined grained view\n"))
#********************

quant.1<-quantile(merge_df_blastn2.uni.fil.1$diff,probs=seq(0,1,0.005),na.rm=TRUE)

#********************
cat(paste(Sys.time()," : classifying rows falling in quantiles >95.5 % or <95.5 %\n"))
#********************

merge_df_blastn2.uni.fil.1$quant<-rep(1,nrow(merge_df_blastn2.uni.fil.1))
merge_df_blastn2.uni.fil.1$quant[merge_df_blastn2.uni.fil.1$diff>quant.1[192]]<-2

#********************
cat(paste(Sys.time()," : adding a column which tells us whether the product annotations are the same or not\n"))
#********************

merge_df_blastn2.uni.fil.1$product.same<-(merge_df_blastn2.uni.fil.1$product.y==merge_df_blastn2.uni.fil.1$product.x)

#********************
cat(paste(Sys.time()," : Do proteins whose genes not localised in the LR genome tend to have different annotations?\n"))
cat(paste(Sys.time()," : % of dissimilar annotation for localised gene pairs: ",(table(merge_df_blastn2.uni.fil.1$product.same,merge_df_blastn2.uni.fil.1$quant)[1]/table(merge_df_blastn2.uni.fil.1$product.same,merge_df_blastn2.uni.fil.1$quant)[2])*100, ", % of dissimilar annotation for distal gene pairs: ",(table(merge_df_blastn2.uni.fil.1$product.same,merge_df_blastn2.uni.fil.1$quant)[3]/table(merge_df_blastn2.uni.fil.1$product.same,merge_df_blastn2.uni.fil.1$quant)[4])*100, "\n"))
cat(paste(Sys.time()," : GENELENGTH RATIO BREAKDOWN?\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio >0.99 and <1.01: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>0.99 & plr <1.01))[1],"\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio >0.95 and <1.05: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>0.95 & plr <1.05))[1],"\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.95 or >=1.05: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.95 | plr>=1.05)),"\n"))
cat(paste(Sys.time()," : number of localised protein pairs with genelength ratio >0.99 and <1.01: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>0.99 & plr <1.01 & quant==1))[1],"\n"))
cat(paste(Sys.time()," : number of genes in Lr-chr corresponding to the above: "))
str(unique(as.character(subset(merge_df_blastn2.uni.fil.1, plr>0.99 & plr <1.01 & quant==1)[,3])))
cat("\n")
cat(paste(Sys.time()," : number of similar and dissimilar annotations for the above: ",table(subset(merge_df_blastn2.uni.fil.1, plr>0.99 & plr <1.01 & quant==1)$product.same)[2],"and",table(subset(merge_df_blastn2.uni.fil.1, plr>0.99 & plr <1.01 & quant==1)$product.same)[1],"\n"))

cat(paste(Sys.time()," : number of distal protein pairs with genelength ratio >0.95 and <1.05: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>0.95 & plr <1.05 & quant==1))[1],"\n"))
cat(paste(Sys.time()," : number of genes in Lr-chr corresponding to the above: "))
str(unique(as.character(subset(merge_df_blastn2.uni.fil.1, plr>0.95 & plr <1.05 & quant==1)[,3])))
cat("\n")
cat(paste(Sys.time()," : number of similar and dissimilar annotations for the above: ",table(subset(merge_df_blastn2.uni.fil.1, plr>0.95 & plr <1.05 & quant==1)$product.same)[2],"and",table(subset(merge_df_blastn2.uni.fil.1, plr>0.95 & plr <1.05 & quant==1)$product.same)[1],"\n"))

cat(paste(Sys.time()," : Maximum protein length ratio observed: ", merge_df_blastn2.uni.fil.1$plr[which.max(merge_df_blastn2.uni.fil.1$plr)],"\n"))
cat(paste(Sys.time()," : GENELENGTH RATIO OBSESRVATION COUNTS & CORRESPONDING QUANTILE RANGE??\n"))

cat(paste(Sys.time()," : number of protein pairs with genelength ratio >=1.05: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>=1.05))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio >=2: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>=2))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio >=3: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>=3))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio >=6: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr>=6))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr>=1.05)$quant)[2],"distal protein pairs\n"))


cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.95: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.95))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.95)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.95)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.8: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.8))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.8)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.8)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.75: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.75))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.75)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.75)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.50: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.50))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.50)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.50)$quant)[2],"distal protein pairs\n"))
cat(paste(Sys.time()," : number of protein pairs with genelength ratio <=0.25: ",nrow(subset(merge_df_blastn2.uni.fil.1, plr<=0.25))[1],"\t","including",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.25)$quant)[1],"localised protein pairs\t and",table(subset(merge_df_blastn2.uni.fil.1, plr<=0.25)$quant)[2],"distal protein pairs\n"))
#********************
#return(merge_df_blastn2.uni.fil.1)
GeneExamination(merge_df_blastn2.uni.fil.1,merge_df_blastn2.uni.fil.1.watsonratio.eq1,merge_df_blastn2.uni.fil.1.watsonratio.eq2,merge_df_blastn2.uni.fil.1.watsonratio.gteq2,merge_df_blastn2.uni.fil.1.glginproduct.x,merge_df_blastn2.uni.fil.1.subset_dup.sracfaa,tmp)
}


#igraph
GeneExamination<-function(genedf,genedf.watsonratio.eq1,genedf.watsonratio.eq2,genedf.watsonratio.gteq2,genedf.glginproduct.x,genedf.subset_dup.sracfaa,tmp)
{
#********************
cat(paste(Sys.time()," : EXAMINATION OF CO-LOCALISED GENES IN DETAIL, SEARCHING FOR SPLIT GENES USING iGRAPH\n"))
cat(paste(Sys.time()," : converting genedf to a graph and decomposing into connected components..\n"))
#********************

net<-graph.edgelist(as.matrix(genedf[c(3,2)]))
E(net)$plr<-genedf$plr
E(net)$diff<-genedf$diff
E(net)$product.same<-genedf$product.same
comp.net<-decompose(net)
#********************
cat(paste(Sys.time()," : make a column to capture the connected component memberships\n"))
#********************

comp.net.num.both<-sapply(comp.net,vcount)
comp.net.num.lrac<-sapply(comp.net,FUN=function(x,y){length(intersect(V(x)$name,y))},y=as.character(genedf$lrac_ORFid))
comp.net.num.srac<-sapply(comp.net,FUN=function(x,y){length(intersect(V(x)$name,y))},y=as.character(genedf$bin1_srac_ORFid))
comp.net.summary<-data.frame(ccid=1:length(comp.net),num.nodes=comp.net.num.both,num.lrac=comp.net.num.lrac,num.srac=comp.net.num.srac,stringsAsFactors=F)
#********************
cat(paste(Sys.time()," : extract midpoint difference values..\n"))
#********************

comp.net.diff<-lapply(comp.net,FUN=function(x){E(x)$diff})

#********************
cat(paste(Sys.time()," : calculate the minimum difference within each cluster\n"))
#********************

comp.net.diff.min<-sapply(comp.net.diff,FUN=min)
comp.net.diff.mean<-sapply(comp.net.diff,FUN=mean)
comp.net.diff.max<-sapply(comp.net.diff,FUN=max)
comp.net.diff.nna<-sapply(comp.net.diff,FUN=function(x){sum(is.na(x))})

#********************
cat(paste(Sys.time()," : calculate statistics of plr..\n"))
#********************
comp.net.plr<-lapply(comp.net,FUN=function(x){E(x)$plr})
comp.net.plr.min<-sapply(comp.net.plr,FUN=min)
comp.net.plr.mean<-sapply(comp.net.plr,FUN=mean)
comp.net.plr.max<-sapply(comp.net.plr,FUN=max)

#********************
cat(paste(Sys.time()," : Check plr value corresponding the min diff..\n"))
#********************

comp.net.plr.for.min.diff<-func1(comp.net)
dim(comp.net.plr.for.min.diff)

#--put this all together into an updated version of 'comp.net.summary'
comp.net.summary.1<-data.frame(comp.net.summary,diff.min=comp.net.diff.min,diff.mean=comp.net.diff.mean,diff.max=comp.net.diff.max,num.na=comp.net.diff.nna,plr.min=comp.net.plr.min,plr.mean=comp.net.plr.mean,plr.max=comp.net.plr.max,plr.from.min.diff=comp.net.plr.for.min.diff,stringsAsFactors=F)

#--and export

write.table(comp.net.summary.1,file="comp.net.summary.1.txt",sep="\t",col.names=T,row.names=F)
source(file="Plot.r")
GeneQuality_plots(genedf,genedf.watsonratio.eq1,genedf.watsonratio.eq2,genedf.watsonratio.gteq2,genedf.glginproduct.x,genedf.subset_dup.sracfaa,tmp)

}
