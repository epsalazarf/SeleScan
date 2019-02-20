# Selection Scan Pipeline: SFS Empirical vs Simulated
# By: Pavel Salazar

#Requires: concat_sfs_all.sh output file
#Requires: concat_sfs_allsim.sh output file

#<INPUT>
# Data files location
empfile <- "SFS.results.adna.ALL.txt"
simfile <- "SFS.sims.ALL.txt"

# Plot parameters
ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
popcol <- c("red","green","blue","magenta")
plotfile <- "SFS.EMPxSIM.abs.png"

#<START>

message("STARTED> SFS: Empirical vs Simulated")

# Libraries
library(data.table)

# Data
message("Reading empirical data...")
message("File: ", empfile)
## Read Empirical Data
empdata <- data.frame(fread(empfile,header=T, drop= c("CHR","POS")))
colnames(empdata) <- c(pops)
empdata <- round(empdata,3)
emptotal <- dim(empdata)[1]
message("Loaded data: ", emptotal, " SNPs found")

## Read simulated data
message("Reading simulated data...")
message("File: ", simfile)
simdata <- data.frame(fread(simfile,header=T, drop= c("SIM","POS")))
colnames(simdata) <- c(pops)
simdata <- round(simdata,3)
simtotal <- dim(simdata)[1]
message("Loaded data: ", simtotal, " SNPs found")
if(max(simdata) > 0.5){
  message("MAFs >0.5 found: folding frequencies...")
  simdata <- data.frame(lapply(simdata, function(x) ifelse(x>0.5,abs(x-1),x)))
}

#<CORE>

#Remove singletons from simdata
message("Removing singletons from simulated data...")
singls <- sapply(lapply(apply(simdata,2, unique),sort),"[[",2)
szrows <- c()
for(i in 1:length(pops)){szrows <- c(szrows,intersect(which(simdata[,i] == singls[i]),which(apply(simdata[,-i],1,sum) == 0)))}
simdata <- simdata[-szrows,]
simtotal <- dim(simdata)[1]
message("Singletons filtered: ", simtotal, " SNPs remain")

# True Binning
message("Creating bins for data...")
empbins <- unique(sort(unlist(apply(empdata,2, unique))))
counte <- data.frame(bin= empbins)
for (p in pops) {
  x <- table(empdata[,p])
  for (i in names(x)) {
    counte[counte$bin == i,p] <- x[i]
  }
}
message("Unique bins for empirical data: ", length(empbins))

simbins <- unique(round(sort(unlist(apply(simdata,2, unique))),3))
counts <- data.frame(bin= simbins)
for (p in pops) {
  x <- table(simdata[,p])
  for (i in names(x)) {
    counts[counts$bin == i,p] <- x[i]
  }
}
message("Unique bins for simulated data: ", length(simbins))

if (all(empbins == simbins)){
  message("Checkpoint: Both datasets have identical bins")
} else {
  message("Warning: differing bins between datasets")
}

#Plot abs
message("Plotting Site Frequency Spectrum...")
png(plotfile, width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
ymin <- floor(min(min(na.omit(counte[,-1])),min(na.omit(counts[,-1])))/10000)*10000
ymax <- ceiling(max(max(na.omit(counte[,-1])),max(na.omit(counts[,-1])))/10000)*10000
plot(c(0.0,0.5),c(ymin,ymax),type="n", xaxt="n", yaxt="n", main= "SFS: Empirical vs Simulated", xlab= "MAF", ylab= "Counts (x1000)")
axis(1, at=seq(0,0.5,0.025),cex.axis=0.8)
axis(2, at=seq(ymin,ymax,10000),cex.axis=0.8,labels= seq(ymin,ymax,10000)/1000,las=1)
for(i in 1:length(pops)){
  lines(counte[!is.na(counte[,i+1]),"bin"],counte[!is.na(counte[,i+1]),pops[i]],col=popcol[i],type="b",pch=16,lwd=2)
}
for(i in 1:length(pops)){
  lines(counts[!is.na(counts[,i+1]),"bin"],counts[!is.na(counts[,i+1]),pops[i]],col=popcol[i],type="b",lwd=1, lty= 3)
}
legend("topright", pops, lty=c(1,1), lwd=c(2,2),col=popcol, cex = 0.8, bty = "n")
dev.off()

if(file.exists(plotfile)){
  message("Success! Plot created")
  message("File:", plotfile)
} else {
  message("Error: Plot not created. Check input.")
}



#Ascertainment Bias Estimation
message("Estimating ascertainment bias...")
ascbtab <- data.frame(bin = empbins)
for (p in pops){
  ascbtab[,p] <- (1 - (counte[,p] / counts[,p]))
}

message("Plotting ascertainment bias...")
png("estimated.ascbias.png", width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
ylims <- c(floor(min(na.omit(ascbtab[,-1]))*10)/10,ceiling(max(na.omit(ascbtab[,-1])*10))/10)
plot(c(0.0,0.5),ylims,type="n", xaxt="n", yaxt="n", main= "Estimated Ascertainment Bias", xlab= "Bin", ylab= "% SNPs ommitted")
axis(1, at=ascbtab$bin,cex.axis=0.7,las=2)
axis(2,at= seq(ylims[1],ylims[2],0.05),cex.axis=0.7,las=1)
for(i in 1:length(pops)){
  lines(ascbtab[!is.na(ascbtab[,i+1]),"bin"],ascbtab[!is.na(ascbtab[,i+1]),pops[i]],col=popcol[i],type="b",pch=16,lwd=2)
}
legend("topright", pops, lty=c(1,1), lwd=c(2,2),col=popcol, cex = 0.8, bty = "n")
dev.off()
if(file.exists("estimated.ascbias.png")){
  message("Success! Plot created")
  message("File: estimated.ascbias.png")
} else {
  message("Error: Plot not created. Check input.")
}

message("ENDED> SFS: Empirical vs Simulated")
#<END>

#<SANDBOX>
# Fixed Percentage Binning
if(0){
  bins <- seq(from = 0, to= 0.5, by= 0.025)
  freqe <- data.frame(bin= bins)
  for (p in pops) {
    x <- table(empdata[,p])/ emptotal
    for (i in names(x)) {
      freqe[max(which(as.numeric(i) >= freqe$bin)),p] <- x[i]
    }
  }
  freqe[is.na(freqe)] <- 0
  rownames(freqe) <- as.character(freqe$bin)
  freqe$bin <- NULL

  freqs <- data.frame(bin= bins)
  for (p in pops) {
    x <- table(simdata[,p])/ simtotal
    for (i in names(x)) {
      freqs[max(which(as.numeric(i) >= freqs$bin)),p] <- x[i]
    }
  }
  freqs[is.na(freqs)] <- 0
  rownames(freqs) <- as.character(freqs$bin)
  freqs$bin <- NULL

  #Plot freqs
  ymax <- round(max(max(na.omit(freqe[,-1])),max(na.omit(freqs[,-1])))*10)/10

  png("ALL.SFS.mixed.lines.png", width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
  plot(c(0.0,0.5),c(0.0,ymax),type="n", xaxt="n", yaxt="n", main= "ALL - SFS (mixed)", xlab= "MAF", ylab= "Freq")
  axis(1, at=seq(0,0.5,0.025),cex.axis=0.8)
  axis(2, at=seq(0,ymax,0.025),cex.axis=0.8)
  for(i in 1:length(pops)){
    lines(freqe[!is.na(freqe[,i+1]),"bin"],freqe[!is.na(freqe[,i+1]),pops[i]],col=popcol[i],type="b",pch=16,lwd=2)
  }
  for(i in 1:length(pops)){
    lines(freqs[!is.na(freqs[,i+1]),"bin"],freqs[!is.na(freqs[,i+1]),pops[i]],col=popcol[i],type="b",lwd=1, lty= 3)
  }
  legend("topright", pops, lty=c(1,1), lwd=c(2,2),col=popcol, cex = 0.8)
  dev.off()
}
