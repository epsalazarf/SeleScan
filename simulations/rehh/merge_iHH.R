#Activated by RUN_PIPELINE_2.sh

library("rehh", quietly=T)

args <- commandArgs(trailingOnly = TRUE)

file<-args[1]
replicates_neutral<-args[2]

if(file.exists(paste("ihh/sim_1.",file,".ihh",sep=""))){
	load(file=paste("ihh/sim_1.",file,".ihh",sep=""))
}else{
	load(file=paste("ihh/sim_2.",file,".ihh",sep=""))
}

ihh_neutral<-ihh_neutral_temp

for (i in 2:replicates_neutral) {
	print(i)
	if(file.exists(paste("ihh/sim_",i,".",file,".ihh",sep=""))){
		load(file=paste("ihh/sim_",i,".",file,".ihh",sep=""))
		ihh_neutral<-rbind(ihh_neutral,ihh_neutral_temp)
	}
}

save(ihh_neutral,file=paste("ihh/",file,".ihh",sep=""))
