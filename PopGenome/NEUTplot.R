# Selection Scan Pipeline: Neutrality Plots
# By: Pavel Salazar

#Requires: concat_neut_per_pop.sh output files


#<INPUT>

ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
prefix <-("results/neut.results.adna.")
stat <- "TD" #TD: Tajima's D, FLF,FLD: Fu and Li's F & D

#<START>
message("STARTED> Neutrality Plot")

# Libraries
library(data.table)

# UCSC Chromosome Palette: (Luminosity 80%)
chrcols <- c("#f4bc7f","#caca88","#cbcb67","#ff8989","#ff7474","#ff77ff","#edbaba","#ffae3e","#f3c000","#cece00","#a6d800","#aef49b","#92d87e","#babaff","#98c4ff","#98cbfe","#00e0e0","#9cd0d0","#faa1ff","#ff92ff","#e3b1ff","#c6c6c6")

#<CORE>
for (p in pops) {

  #Load and clean data
  message("Reading data for ", p,"...")
  file <- (paste(prefix,p,".txt", sep= ""))
  results <- data.frame(fread(file, header=T))
  results <- results[complete.cases(results),]

  #Top Tables
  stat99 <- c(quantile(results[,stat],prob=0.005, na.rm = T),quantile(results[,stat],prob=0.995, na.rm = T))
  stat999 <- c(quantile(results[,stat],prob=0.0005, na.rm = T),quantile(results[,stat],prob=0.9995, na.rm = T))
  restop1 <- subset(results, results[,stat] > stat99[2] | results[,stat] < stat99[1] )
  write.table(restop1,file = paste("results/",p,".",stat,".top1.txt",sep = ""), quote = F, sep = "\t", row.names = F)
  message("Wrote Top 1% ",  stat, " table for ", p)
  restop01 <- subset(results, results[,stat] > stat999[2] | results[,stat] < stat999[1])
  write.table(restop01,file = paste("results/",p,".",stat,".top01.txt",sep = ""), quote = F, sep = "\t", row.names = F)
  message("Wrote Top 0.1% ",  stat, " table for ", p)

  # Plotting
  message("Plotting data for ", p,"...")
  results$Col <- chrcols[results$CHR]

  # Generate chromosome labels
  results$genpos=NA
  lastbase=0
  ticks=NULL
  for (i in unique(results$CHR)) {
    if (i==1) {
      results[results$CHR==i, ]$genpos=results[results$CHR==i, ]$start
    } else {
      lastbase=lastbase+tail(subset(results,CHR==i-1)$start, 1)
      results[results$CHR==i, ]$genpos=results[results$CHR==i, ]$start+lastbase
    }
    ticks = c(ticks, (min(results[results$CHR == i,]$genpos) + max(results[results$CHR == i,]$genpos))/2 + 1)
  }
  labs <- unique(results$CHR)

  # Plot parameters
  xmax = ceiling(max(results$genpos) * 1.03)
  xmin = floor(max(results$genpos) * -0.03)
  ymin <- floor(min(results[,stat]))
  ymax <- ceiling(max(results[,stat]))
  def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="Chromosome", ylab=stat, cex.axis=0.7,cex = 0.7,main= p)

  # Generate plot
  png(filename = paste("results/",p,".",stat,".plot.png",sep=""), width = 1920, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
  do.call("plot", c(NA, def_args))
  axis(1, at=ticks, labels=labs)
  axis(2, at=seq(ymin,ymax,0.5),cex.axis=0.8)

  # Plot points (SLOW)
  icol=1
  for (i in unique(results$CHR)) {
    #with(results[results$CHR==i,], points(genpos, TD, col=chrcols[icol], pch=20))
    temp <- subset(results,results$CHR==i)
    points(temp$genpos, temp[,stat], col=chrcols[icol], pch=20)
    icol=icol+1
  }
  abline(h = 0, col ="grey" , lwd = 2)
  abline(h = stat99, col = "blue", lwd = 1)
  abline(h = stat999, col = "red", lwd = 1)
  legend("topright", c("99%","99.9%"), lty=c(1,1), lwd=c(1,1),col=c("blue", "red"), cex = 0.6, bty = "n")
  dev.off()
  message(stat," plotted for ", p)
  message("Finished: ",p)
}

message("Finished all tasks!")

message("ENDED> Neutrality Plot")
#<END>
