#Activated by run_conversion_and_parsing.pl

library("rehh", quietly=T)

####

args <- commandArgs(trailingOnly = TRUE)

file<-args[1]

print(file)

hap<-data2haplohh(hap_file=paste("input/",file,".hap",sep=""),map_file=paste("input/",file,".map",sep=""))

ihh_neutral_temp<-scan_hh(hap)

save(ihh_neutral_temp,file=paste("ihh/",file,".ihh",sep=""))
