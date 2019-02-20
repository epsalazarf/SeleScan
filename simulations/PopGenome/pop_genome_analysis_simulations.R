#Activated by RUN_PopGenome.sh

# <START>

# <DATA>

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
seq_len <- args[2]

message("Starting neutrality tests for ",input_file)

library(methods, quietly=T)
library("PopGenome", quietly=T)

sim.data <- readMS(file=paste("../input_data/",input_file,".ms",sep=""))

positions <- read.table(paste("../input_data/",input_file,".ms.pos",sep=""),header=T)
positions <- round(positions[,1] * as.integer(seq_len))

sim.data@region.data@biallelic.sites[[1]]<-positions

sim.data@n.sites[[1]]<-as.integer(seq_len)

# </DATA>

# <OUTGROUP>

outgroup <- c("9")
sim.data <- set.outgroup(sim.data,new.outgroup=outgroup,diploid=FALSE)

#ALL
POP1 <- as.character(c(1:2))

POP2 <- as.character(c(3:4))

POP3 <- as.character(c(5:6))

POP4 <- as.character(c(7:8))

sim.data <- set.populations(sim.data,list(POP1,POP2,POP3,POP4), diploid=F)

# <SLIDING.WINDOWS>

sim.data.slide <- sliding.window.transform(sim.data,width=25000, jump=12500,type=2,whole.data=TRUE)

# <NEUTRALITY.TESTS>


sim.data.slide <- neutrality.stats(sim.data.slide,FAST=TRUE)
neutral_tests.sim <- get.neutrality(sim.data.slide)


message("Finished neutrality tests for ",input_file)

message("Generating output...")

write.table(neutral_tests.sim[1],file=paste("results/neutrality.",input_file,".POP1.txt",sep=""))
write.table(neutral_tests.sim[2],file=paste("results/neutrality.",input_file,".POP2.txt",sep=""))
write.table(neutral_tests.sim[3],file=paste("results/neutrality.",input_file,".POP3.txt",sep=""))
write.table(neutral_tests.sim[4],file=paste("results/neutrality.",input_file,".POP4.txt",sep=""))

message("Starting Fst and SFS for ",input_file)

#  <MAF> (SNP-based)
sim.data <- detail.stats(sim.data, site.FST=TRUE)
write.table(t(sim.data@region.stats@minor.allele.freqs[[1]]),file=paste("results/SFS.",input_file,".ALL.txt",sep=""))

# <FST> (SNP-based)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".ALL.txt",sep=""))

# MAF per POP
sim.data <- detail.stats(sim.data,new.populations=list(POP1))
write.table(t(sim.data@region.stats@minor.allele.freqs[[1]]),file=paste("results/SFS.",input_file,".POP1.txt",sep=""))

sim.data <- detail.stats(sim.data,new.populations=list(POP2))
write.table(t(sim.data@region.stats@minor.allele.freqs[[1]]),file=paste("results/SFS.",input_file,".POP2.txt",sep=""))

sim.data <- detail.stats(sim.data,new.populations=list(POP3))
write.table(t(sim.data@region.stats@minor.allele.freqs[[1]]),file=paste("results/SFS.",input_file,".POP3.txt",sep=""))

sim.data <- detail.stats(sim.data,new.populations=list(POP4))
write.table(t(sim.data@region.stats@minor.allele.freqs[[1]]),file=paste("results/SFS.",input_file,".POP4.txt",sep=""))

# Paired FST

sim.data<-detail.stats(sim.data, new.populations=list(POP1,POP2), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP1xPOP2.txt",sep=""))

sim.data<-detail.stats(sim.data, new.populations=list(POP1,POP3), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP1xPOP3.txt",sep=""))

sim.data<-detail.stats(sim.data, new.populations=list(POP1,POP4), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP1xMAY.txt",sep=""))

sim.data<-detail.stats(sim.data, new.populations=list(POP2,POP3), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP2xPOP3.txt",sep=""))

sim.data<-detail.stats(sim.data, new.populations=list(POP2,POP4), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP2xPOP4.txt",sep=""))

sim.data<-detail.stats(sim.data, new.populations=list(POP3,POP4), site.FST=TRUE)
write.table(sim.data@region.stats@site.FST[[1]],file=paste("results/fst.results.",input_file,".POP3xPOP4.txt",sep=""))

message("Finished Fst and SFS for ",input_file)

# <END>
