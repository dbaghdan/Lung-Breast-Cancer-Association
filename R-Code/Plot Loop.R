#Creates a Manhattan and QQ plot for the association file of each tissue on a list

#Open Relevant Libraries
library(dplyr)
library(ggplot2)
library(sas7bdat)
library(tidyr)
library(qqman)
#Read in the Chromosome, Ensemble, and gene names and change the column names
Chrome=read.table("~/Desktop/ChrENGene.txt")
colnames(Chrome)<-c('CHR','gene','gene_name')
#Change the chromosomes into numbers
Chrome<-transform(Chrome, CHR=as.numeric(CHR))
#Read in the Tissue list
x<-scan("~/Desktop/predixcan/dbs.txt",what="",sep="\n")
x<-strsplit(x,"[[:space:]]+")
#For each tisssue in the list
for(item in x){
#Read in the Data
  s1="~/Desktop/predixcan/"
  s2=item
  s3="/association_0.1.txt"
  path=paste(s1,s2,s3,sep="")
  Tiss=read.table(path,header=T)
  Tissu=na.omit(Tiss)
  #Add the Chromosome and gene names
  Tissue<-left_join(Tissu,Chrome, by='gene')
  s1="~/Desktop/PredixcanPlots/"
  s3="_qqplot_0.1.jpg"
  path=paste(s1,s2,s3,sep="")
  #Open access to the desire jpeg location
  jpeg(file = path)
  #Write a qq plot
  qq(Tissue$p, main='QQplot of Tissue', psh=18)
  #Close
  dev.off()
  s3="_manplot_0.1.jpg"
  path=paste(s1,s2,s3,sep="")
  #Open desired jpeg location
  jpeg(filename = path)
  #Write manhattan plot, labelling genes at a predetermined threshold 
  print(ggplot(Tissue, aes(x=CHR, y=-log10(p),color=CHR)) +geom_point(aes(size = 10,)) + geom_text(aes(label=ifelse(-log10(p)>4.7,as.character(gene_name),'')),hjust= 1.2, vjust=0) + ggtitle("Manhattan Plot Tissue") + xlab("Chromosome") + ylab("-log10(p)") + theme_set(theme_gray(base_size = 18))+ theme(legend.position="none") + scale_x_continuous(breaks=seq(1,25,1)))
 #Close
  dev.off()
  }
