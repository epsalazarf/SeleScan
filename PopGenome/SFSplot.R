# Selection Scan Pipeline: SFS Plot
# By: Pavel Salazar

#Requires: concat_sfs_per_pop.sh output files

#<INPUT>

ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
prefix <-("results/neut.results.adna.")
popcol <- c("red","green","blue","magenta")
plots <- c("BB","UL") # BB: binned bars, UL: unbinned lines

#<START>
message("STARTED> SFS Plot")

# Libraries
library(data.table)

# Data: merge multiple SFS tables
results <- data.frame()
for (p in pops) {
  message("Reading data for ", p,"...")
  if(p==pops[1]){
    results <- data.frame(fread(paste(prefix,p,".txt", sep= ""),header=T))
  } else {
    results <- cbind(results,data.frame(fread(paste(prefix,p,".txt", sep= ""),header=T, drop= c("CHR","POS"))))
  }
}
colnames(results) <- c("CHR","POS",pops)
total <- dim(results)[1]

#<CORE>

# BB: Binned Barplot
if("BB" %in% plots){
  freqs <- data.frame(bin= seq(from = 0, to= 0.5, by= 0.025))
  for (p in pops) {
    x <- table(results[,p])/ total
    for (i in names(x)) {
      freqs[max(which(as.numeric(i)+0.01 >= freqs[,1])),p] <- x[i]
    }
  }
  freqs[is.na(freqs)] <- 0
  rownames(freqs) <- as.character(freqs$bin)
  freqs$bin <- NULL
  png("results/ALL.SFS.bars.png", width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
  barplot(t(freqs),beside=T, main= "ALL - SFS", xlab= "MAF <= X", ylab= "Freq", col=popcol, ylim = c(0,0.2), cex.axis = 0.8, cex.names = 0.8)
  legend("topright", pops, cex=0.8, fill=popcol)
  dev.off()
  message("Output: Binned Barplot")
}


# UL: Unbinned Lineplot
if("UL" %in% plots){
  for (p in pops){
    temp <- data.frame(table(results[,p])/total)
    colnames(temp) <- c("freq",p)
    if(p==pops[1]){
      freqp <- data.frame(bin= temp[,1])
    }
    freqp <- merge(freqp,temp,by.x=c("bin"), by.y=c("freq"), all= T)
  }

  freqp$bin <- as.numeric(levels(freqp$bin))
  freqp <- freqp[order(freqp$bin),]
  rownames(freqp) <- as.character(round(freqp$bin,3))

  png("results/ALL.SFS.lines.png", width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
  plot(c(0.0,0.5),c(0.0,0.2),type="n", xaxt="n", yaxt="n", main= "ALL - SFS", xlab= "MAF", ylab= "Freq")
  axis(1, at=seq(0,0.5,0.025),cex.axis=0.8)
  axis(2, at=seq(0,0.2,0.025),cex.axis=0.8)
  for(i in 1:length(pops)){
  lines(freqp[!is.na(freqp[,i+1]),"bin"],freqp[!is.na(freqp[,i+1]),pops[i]],col=popcol[i],lwd=2.5)
  points(freqp[!is.na(freqp[,i+1]),"bin"],freqp[!is.na(freqp[,i+1]),pops[i]],col=popcol[i],pch=16)
  }
  legend("topright", pops, lty=c(1,1), lwd=c(2.5,2.5),col=popcol, cex = 0.8)
  dev.off()
  message("Output: Unbinned Lineplot")
}

message("Finished all tasks!")

message("ENDED> SFS Plot")
#<END>
