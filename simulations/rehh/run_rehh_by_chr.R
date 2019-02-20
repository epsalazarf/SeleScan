library("rehh", quietly=T)

args <- commandArgs(trailingOnly = TRUE)

file<-args[1]

hap<-data2haplohh(hap_file=paste("input/",file,".hap",sep=""),map_file=paste("input/",file,".map",sep=""))

pdf(paste("plots/",file,".pdf",sep=""))

res<-scan_hh(hap)
ihs<-ihh2ihs(res)
distribplot(ihs$res.ihs[,3])
ihsplot(ihs$res.ihs)
dev.off()

write.table(ihs$res.ihs,file=paste("results/",file,".ihs_temp",sep=""))

print(warnings())
