#Goes through each tissue's association file and filters based on R^2 value

#Open relevant Libraries
library(dplyr)
library(sas7bdat)
library(tidyr)
#Read in the list of tissues and run a loop
x<-scan("~/Desktop/predixcan/dbsr2.txt",what="",sep="\n")
x<-strsplit(x,"[[:space:]]+")
#For each tissue in the list
for(item in x){
#Read in the data table
  s1="~/Desktop/predixcan/"
  s2=item
  s3="0_0.5.db/association.txt"
  path=paste(s1,s2,s3,sep="")
  Tiss=read.table(path,header=T)
  Tiss=na.omit(Tiss)
  #Read in and append the R^2 values to the table
  s1="~/mount/wheelerlab1/Data/PrediXcan_R2/"
  s3="0.5_hapmapSnpsCEU_unscaled.allResults.txt"
  path=paste(s1,s2,s3,sep="")
  Rvalues=read.table(path,header=T)
  FiltR<-left_join(Tiss,Rvalues, by='gene')
  #Filter the genes by a specified R^2 value
  FiltR<-dplyr::filter(FiltR,R2>=0.1)
  #Write the new file
  s1="~/Desktop/predixcan/"
  s3="0_0.5.db/association_0.1.txt"
  path=paste(s1,s2,s3,sep="")
  write.table(FiltR,path,quote=FALSE,col.names=TRUE)
}
