#Goes through the association file of each tissue on a list and filters based on P value

#Open relevant libraries
library(dplyr)
library(sas7bdat)
library(tidyr)
#Run the program fto read in the first tissue
path="~/Desktop/predixcan/DGN-WB-unscaled_0.5.db/association.txt"
Chrome=read.table(path,header=T)
#Remove NA's and add a column for tissue type
Chrome=na.omit(Chrome)
Chrome$Tissue <- rep("DGN",nrow(Chrome))
#Read in a list of all the other tissues and run the loop through each one
x<-scan("~/Desktop/predixcan/dbs.txt",what="",sep="\n")
x<-strsplit(x,"[[:space:]]+")
for(item in x){
  s1="~/Desktop/predixcan/"
  s2=item
  s3="/association.txt"
  path=paste(s1,s2,s3,sep="")
  Tiss=read.table(path,header=T)
  Tiss=na.omit(Tiss)
  Tiss$Tissue <- rep(item,nrow(Tiss))
  #Append each tissue to a compilation table
  Chrome<-rbind(Chrome,Tiss)
}
#Filter the table to get a list of genes with a specified p value
Chrome<-dplyr::filter(Chrome,p<=1e-5)
