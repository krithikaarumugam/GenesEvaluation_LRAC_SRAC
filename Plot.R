

#draft - still refining..
GeneQuality_plots<-function(GQdf,GQdf.watsonratio.eq1,GQdf.watsonratio.eq2,GQdf.watsonratio.gteq2,GQdf.glginproduct.x,GQdf.subset_dup.sracfaa,tmp)
{
plot(GQdf$mplracfaa,GQdf$mpsracfaa.strandnew,col = ifelse(GQdf$sstart.y > GQdf$send.y ,'orange','blue'),pch=16,cex=0.5)

#print("hello")


pdf("plots/mainplot.pdf",8,8)
plot(GQdf$mplracfaa,GQdf$mpsracfaa.strandnew,col = ifelse(GQdf$sstart.y > GQdf$send.y ,'orange','blue'),pch=16,cex=0.5)
dev.off()

#print("hello")
print(head(GQdf$mplracfaa))

#plot difference between midpoints Vs genelength ratio
plot((GQdf$mplracfaa-GQdf$mpsracfaa.strandnew),(GQdf$lrac_ORFlen/GQdf$bin_srac_ORFlen))
#print("hello")


#check diamond stats of the 2271 genes
plot((GQdf.watsonratio.eq1$lrac_ORFlen/GQdf.watsonratio.eq1$bin_srac_ORFlen),GQdf.watsonratio.eq1$pident.x)
text((GQdf.watsonratio.eq1$lrac_ORFlen/GQdf.watsonratio.eq1$bin_srac_ORFlen),GQdf.watsonratio.eq1$pident.x,labels=GQdf.watsonratio.eq1$product.y,cex=0.5,pos=4)
hist(GQdf.watsonratio.eq1$pident.x)
#print("hello")

#mmm, so even though ratio is eq to 1 pid can be less, few hypothetical proteins and one 50S ribosomal protein<90% 
#try plotting midpoints?
plot(GQdf.watsonratio.eq1$mplracfaa,GQdf.watsonratio.eq1$mpsracfaa.strandnew)
#print("hello")

text(GQdf.watsonratio.eq1$mplracfaa,GQdf.watsonratio.eq1$mpsracfaa.strandnew,cex=0.3,labels=GQdf.watsonratio.eq1$product.x,pos=4)
#mmm, off diagonal points are hypothetical mostly? maybe there's another way to plot?
#plot genelength ratio and midpoint difference
plot((GQdf.watsonratio.eq1$mplracfaa-GQdf.watsonratio.eq1$mpsracfaa.strandnew),(GQdf.watsonratio.eq1$lrac_ORFlen/GQdf.watsonratio.eq1$bin_srac_ORFlen))
#print("hello")


#check diamond stats
plot((GQdf.watsonratio.gteq2$lrac_ORFlen/GQdf.watsonratio.gteq2$bin_srac_ORFlen),GQdf.watsonratio.gteq2$pident.x)
text((GQdf.watsonratio.gteq2$lrac_ORFlen/GQdf.watsonratio.gteq2$bin_srac_ORFlen),GQdf.watsonratio.gteq2$pident.x,labels=GQdf.watsonratio.gteq2$product.x,cex=0.5)
#mostly hypothetical?
#plot their midpoints
plot(GQdf.watsonratio.gteq2$mplracfaa,GQdf.watsonratio.gteq2$mpsracfaa.strandnew)
text(GQdf.watsonratio.gteq2$mplracfaa,GQdf.watsonratio.gteq2$mpsracfaa.strandnew,labels=GQdf.watsonratio.gteq2$product.x,cex=0.5,pos=4)

#plot genelength ratio and midpoint difference
plot((GQdf.watsonratio.gteq2$mplracfaa-GQdf.watsonratio.gteq2$mpsracfaa.strandnew),(GQdf.watsonratio.gteq2$lrac_ORFlen/GQdf.watsonratio.gteq2$bin_srac_ORFlen))
text((GQdf.watsonratio.gteq2$mplracfaa-GQdf.watsonratio.gteq2$mpsracfaa.strandnew),(GQdf.watsonratio.gteq2$lrac_ORFlen/GQdf.watsonratio.gteq2$bin_srac_ORFlen),labels=GQdf.watsonratio.gteq2$product.x,cex=0.5,pos=4)

#glg
#plot genelength ratio and midpoint difference
plot((GQdf.glginproduct.x$mplracfaa-GQdf.glginproduct.x$mpsracfaa.strandnew),(GQdf.glginproduct.x$lrac_ORFlen/GQdf.glginproduct.x$bin_srac_ORFlen))
text((GQdf.glginproduct.x$mplracfaa-GQdf.glginproduct.x$mpsracfaa.strandnew),(GQdf.glginproduct.x$lrac_ORFlen/GQdf.glginproduct.x$bin_srac_ORFlen),labels=GQdf.glginproduct.x$product.x,cex=0.5)
#so all good ?

plot(GQdf.glginproduct.x$mplracfaa,GQdf.glginproduct.x$mpsracfaa.strandnew)



#subset for sracfaa duplicate #plot midpoints, genelength ratio and midpoint difference
plot(GQdf.subset_dup.sracfaa$mplracfaa,GQdf.subset_dup.sracfaa$mpsracfaa.strandnew)
plot((GQdf.subset_dup.sracfaa$mplracfaa-GQdf.subset_dup.sracfaa$mpsracfaa.strandnew),(GQdf.subset_dup.sracfaa$lrac_ORFlen/GQdf.subset_dup.sracfaa$bin_srac_ORFlen))
text((GQdf.subset_dup.sracfaa$mplracfaa-GQdf.subset_dup.sracfaa$mpsracfaa.strandnew),(GQdf.subset_dup.sracfaa$lrac_ORFlen/GQdf.subset_dup.sracfaa$bin_srac_ORFlen),labels=GQdf.subset_dup.sracfaa$product.x,cex=0.3)
#print("hello")

}


