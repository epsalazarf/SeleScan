# Selection Scan Pipeline: Paired Fst Plot
# By: Pavel Salazar

#Requires: concat_all_fst.sh output file

#<INPUT>

prefix <-("results/fst.results.adna.")
ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
pairs <- combn(pops, 2)
paired <- c()
for (i in 1:ncol(pairs)){
  paired <- c(paired,paste0(pairs[1,i],"x",pairs[2,i]))
}

#<START>
message("STARTED> Paired Fst Plot")

# Libraries
library(data.table)

#<CORE>
for(p in paired){
  # Data
  results <- data.frame(fread(paste(prefix,p,".txt", sep=""), header=T))

  #Top Tables
  stat99 <- quantile(results$FST,prob=0.99, na.rm = T)
  stat999 <- quantile(results$FST,prob=0.999, na.rm = T)
  restop1 <- subset(results, results$FST > stat99)
  write.table(restop1,file = paste("results/",p,".fst.top1.txt",sep = ""), quote = F, sep = "\t", row.names = F)
  message("Wrote Top 1% Fst table for ", p)
  restop01 <- subset(results, results$FST > stat999)
  write.table(restop01,file = paste("results/",p,".fst.top01.txt",sep = ""), quote = F, sep = "\t", row.names = F)
  message("Wrote Top 0.1% Fst table for ", p)

  # Plot
  message("Plotting data for ", p,"...")

  # UCSC Chromosome Palette: (Luminosity 80%)
  chrcols <- c("#f4bc7f","#caca88","#cbcb67","#ff8989","#ff7474","#ff77ff","#edbaba","#ffae3e","#f3c000","#cece00","#a6d800","#aef49b","#92d87e","#babaff","#98c4ff","#98cbfe","#00e0e0","#9cd0d0","#faa1ff","#ff92ff","#e3b1ff","#c6c6c6")
  results$Col <- chrcols[results$CHR]

  # Generate chromosome labels
  results$genpos=NA
  lastbase=0
  ticks=NULL
  for (i in unique(results$CHR)) {
    if (i==1) {
      results[results$CHR==i, ]$genpos=results[results$CHR==i, ]$POS
    } else {
      lastbase=lastbase+tail(subset(results,CHR==i-1)$POS, 1)
      results[results$CHR==i, ]$genpos=results[results$CHR==i, ]$POS+lastbase
    }
    ticks = c(ticks, (min(results[results$CHR == i,]$genpos) + max(results[results$CHR == i,]$genpos))/2 + 1)
  }
  labs <- unique(results$CHR)

  # Plot parameters
  xmax = ceiling(max(results$genpos) * 1.03)
  xmin = floor(max(results$genpos) * -0.03)
  ymin <- round((min(results$FST)-0.024)/0.05)*0.05
  ymax <- round((max(results$FST)+0.024)/0.05)*0.05
  def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="Chromosome", ylab="FST", cex.axis=0.7,cex = 0.7,main= p)

  # Generate plot
  png(filename = paste("results/",p,".Fst.plot.png",sep=""), width = 1920, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
  do.call("plot", c(NA, def_args))
  axis(1, at=ticks, labels=labs)
  axis(2, at=seq(ymin,ymax,0.05),cex.axis=0.8)
  icol=1
  for (i in unique(results$CHR)) {
    with(results[results$CHR==i,], points(genpos, FST, col=chrcols[icol], pch=20))
    icol=icol+1
  }
  abline(h = 0, col ="grey" , lwd = 2)
  abline(h = stat99, col = "blue", lwd = 1)
  abline(h = stat999, col = "red", lwd = 1)
  legend("topright", c("99%","99.9%"), lty=c(1,1), lwd=c(1,1),col=c("blue", "red"), cex = 0.6, bty = "n")
  dev.off()
  message("Paired Fst plotted for ", p)
}
message("Finished all tasks")

message("ENDED> Fst Plot")
#<END>
