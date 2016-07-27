#Takes two sets of Predixcan Association Data and creates a merged QQ plot

#Open relevant Libraries
library(dplyr)
library(ggplot2)
library(sas7bdat)
library(tidyr)
library(qqman)
#Open qqunif
source("~/mount/wheelerlab1/Data/qq-man-scripts/qqunif.r")
#Read in the Chromosomes, ensembles, and gene names, and change column names
Chrome=read.table("~/Desktop/ChrENGene.txt")
colnames(Chrome)<-c('CHR','gene','gene_name')
#Turn chromosomes into numerics
Chrome<-transform(Chrome, CHR=as.numeric(CHR))
#Read in the list of tissues
x<-scan("~/Desktop/predixcan/DB3.txt",what="",sep="\n")
x<-strsplit(x,"[[:space:]]+")
#For each tissue
for(item in x){
#Read in the Association Data
  s1="~/Desktop/predixcan/"
  s2=item
  s3="/association_0.01.txt"
  path=paste(s1,s2,s3,sep="")
  Tiss=read.table(path,header=T)
  Tissu=na.omit(Tiss)
  #Merge to add Chromosome and gene name
  Tissue<-left_join(Tissu,Chrome, by='gene')
  #Open the association from the second predixcan
  s1="Downloads/Pred_assoc/"
  s3=".txt"
  path=paste(s1,s2,s3,sep="")
  Org=read.table(path,header=T)
  Org=na.omit(Org)
  #Filter for p value
  Org2=filter(Org,p<=0.05)
  #Merge the two scans
  Tissue2=Tissue[Tissue$gene %in% Org2$gene,]
  nn=length(Tissue2$p)
  xx= -log10((1:nn)/(nn+1))
  s1="~/Desktop/PredixcanPlots/MergedPlots/"
  s3="_mergeqq_0.01.jpg"
  path=paste(s1,s2,s3,sep="")
  #Sort the data by p value
  Tissue2=Tissue2[rev(order(Tissue2$p)),]
  #Open the desire path 
  jpeg(file = path)
  #Write a qqplot for one PrediXcan
  qqunif(Tissue$p,plot=T, BH=F)
  #Add the points for the other one
  points(sort(xx), sort(-log10(Tissue2$p)), col="red")
  #Filter by desired p value
  Tissue2=filter(Tissue2,p<=0.005)
  nn=length(Tissue2$p)
  xx= -log10((1:nn)/(nn+1))
  #Label points from the second PrediXcan by gene name (Use adj to move the labels horizontally
  if(nrow(Tissue2)!=0){
    text(sort(xx), y=sort(-log10(Tissue2$p)), labels = (Tissue2$gene), adj = -1.5)
  }
  #Add legend
  legend('bottomright', c("Genes for Breast Cancer", "All Genes"), col=c("red","black"), pch=19,cex=1)
  dev.off()
}
