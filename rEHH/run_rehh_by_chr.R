#Activated by run_conversion_and_parsing.pl

library("rehh", quietly=T)
source("rehh1devplot.R")

args <- commandArgs(trailingOnly = TRUE)

file<-args[1]

hap<-data2haplohh(hap_file=paste(file,".hap",sep=""),map_file=paste(file,".map",sep=""))

pdf(paste(file,".pdf",sep=""))

res<-scan_hh(hap)
ihs<-ihh2ihs(res,freqbin=0.05)
distplotmin(ihs$iHS[,3])
ihsplotmin(ihs)
dev.off()

write.table(ihs$iHS,file=paste(file,".ihs_temp",sep=""))

print(warnings())
