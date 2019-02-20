# Selection Scan Pipeline: XP-EHH Plot
# By: Pavel Salazar

#Requires:  run_conversion_and_parsing.pl output files

#<INPUT>
ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
chrms <- c(1:22)
output <- c("T","P") #T: tables, P: plots

#<START>
message("STARTED> XP-EHH Plot")

# Libraries
library("rehh")
library("qqman")

# Functions
antilog10<-function(lx){exp(lx/log(exp(1),10))}

# Preparations
pairs <- combn(pops, 2)
## UCSC Chromosome Palette: (Luminosity 80%)
chrcols <- c("#f4bc7f","#caca88","#cbcb67","#ff8989","#ff7474","#ff77ff","#edbaba","#ffae3e","#f3c000","#cece00","#a6d800","#aef49b","#92d87e","#babaff","#98c4ff","#98cbfe","#00e0e0","#9cd0d0","#faa1ff","#ff92ff","#e3b1ff","#c6c6c6")
names(chrcols) <- c(1:22)

#<CORE>
for (p in 1:ncol(pairs)){
  pairname <-  paste0(pairs[1,p],"x",pairs[2,p])
  message("Reading data for ",pairname,"...")
  for(i in chrms){
    file <- paste0("chr",i,"/",pairs[1,p],".chr",i)
    hap<-data2haplohh(paste0(file,".hap"),paste0(file,".map"),chr.name=i)
    res<-scan_hh(hap)
    if(i==1){wg.res.pop1<-res}else{wg.res.pop1<-rbind(wg.res.pop1,res)}
    file <- paste0("chr",i,"/",pairs[2,p],".chr",i)
    hap<-data2haplohh(paste0(file,".hap"),paste0(file,".map"),chr.name=i)
    res<-scan_hh(hap)
    if(i==1){wg.res.pop2<-res}else{wg.res.pop2<-rbind(wg.res.pop2,res)}
  }
  message("Intersecting data for ",pairname,"...")
  sharedsnps <- intersect(row.names(wg.res.pop1),row.names(wg.res.pop2))
  wg.res.pop1 <- wg.res.pop1[sharedsnps,]
  wg.res.pop2 <- wg.res.pop2[sharedsnps,]
  wg.xpehh<-ies2xpehh(wg.res.pop1,wg.res.pop2)
  wg.xpehh$ID <- row.names(wg.xpehh)

  #<OUTPUT>

  # Top Tables
  if("T" %in% output){
    write.table(wg.xpehh[,c("ID","CHR","POSITION","XPEHH","-log10(p-value) [bilateral]")],file = paste0("results/",pairname,".xpehh.all.txt"), quote = F, sep = "\t", row.names = F)
    message("Wrote ALL table for ", pairname)
    restop1 <- subset(wg.xpehh,wg.xpehh$XPEHH > quantile(wg.xpehh$XPEHH,prob=0.99, na.rm = T))
    write.table(restop1,file = paste0("results/",pairname,".xpehh.top1.txt"), quote = F, sep = "\t", row.names = F)
    message("Wrote Top 1% table for ", pairname)
    restop01 <- subset(wg.xpehh,wg.xpehh$XPEHH > quantile(wg.xpehh$XPEHH,prob=0.999, na.rm = T))
    write.table(restop01,file = paste0("results/",pairname,".xpehh.top01.txt"), quote = F, sep = "\t", row.names = F)
    message("Wrote Top 0.1% table for ", pairname)
  }

  # Plotting
  if("P" %in% output){
    message("Transforming data for plotting ",pairname,"...")
    wg.xpehh$CHR <- as.numeric(wg.xpehh$CHR)
    wg.xpehh$ID <- row.names(wg.xpehh)
    wg.xpehh$Praw <- antilog10(wg.xpehh[,4])
    wg.xpehh$Col <- chrcols[wg.xpehh$CHR]

    # Generate chromosome labels
    wg.xpehh$genpos=NA
    lastbase=0
    ticks=NULL
    for (j in unique(wg.xpehh$CHR)) {
      if (j==1) {
        wg.xpehh[wg.xpehh$CHR==j, ]$genpos=wg.xpehh[wg.xpehh$CHR==j, ]$POSITION
      } else {
        lastbase=lastbase+tail(subset(wg.xpehh,CHR==j-1)$POSITION, 1)
        wg.xpehh[wg.xpehh$CHR==j, ]$genpos=wg.xpehh[wg.xpehh$CHR==j, ]$POSITION+lastbase
      }
      ticks = c(ticks, (min(wg.xpehh[wg.xpehh$CHR == j,]$genpos) + max(wg.xpehh[wg.xpehh$CHR == j,]$genpos))/2 + 1)
    }
    labs <- unique(wg.xpehh$CHR)
    stat99 <- c(quantile(wg.xpehh$XPEHH,prob=0.005, na.rm = T),quantile(wg.xpehh$XPEHH,prob=0.995, na.rm = T))
    stat999 <- c(quantile(wg.xpehh$XPEHH,prob=0.0005, na.rm = T),quantile(wg.xpehh$XPEHH,prob=0.9995, na.rm = T))

    # Plot parameters
    message("Plotting data for ", pairname,"...")
    xmax = ceiling(max(wg.xpehh$genpos) * 1.03)
    xmin = floor(max(wg.xpehh$genpos) * -0.03)
    ymin <- floor(min(wg.xpehh$XPEHH))
    ymax <- ceiling(max(wg.xpehh$XPEHH))
    def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="Chromosome", ylab="XPEHH", cex.axis=0.7,cex = 0.7,main= paste("XP-EHH:",pairname))

    # Generate plot
    png(filename = paste0("results/",pairname,".xpehh.png"), width = 1920, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
    do.call("plot", c(NA, def_args))
    axis(1, at=ticks, labels=labs)
    axis(2, at=seq(ymin,ymax,0.5),cex.axis=0.8)

    # Plot points (SLOW)
    icol=1
    for (i in unique(wg.xpehh$CHR)) {
      temp <- subset(wg.xpehh,wg.xpehh$CHR==i)
      points(temp$genpos, temp$XPEHH, col=chrcols[icol], pch=20)
      icol=icol+1
    }
    abline(h = 0, col ="grey" , lwd = 2)
    abline(h = stat99, col = "blue", lwd = 1)
    abline(h = stat999, col = "red", lwd = 1)
    legend("topright", c("99%","99.9%"), lty=c(1,1), lwd=c(1,1),col=c("blue", "red"), cex = 0.6, bty = "n")
    dev.off()
    message("XP-EHH scores plotted for ", pairname)

    # qqman: Manhattan Plot for P-values
    pvlines <- c(quantile(wg.xpehh[,4],prob=0.99, na.rm = T),quantile(wg.xpehh[,4],prob=0.999, na.rm = T))
    png(filename = paste0("results/",pairname,".xpehh.logp.png"), width = 1920, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
    manhattan(wg.xpehh, snp = "ID", bp = "POSITION", p= "-log10(p-value) [bilateral]", genomewideline = pvlines[2], suggestiveline = pvlines[1], cex.axis=0.7, cex = 0.7, main= paste("XP-EHH:",pairname), col=chrcols, logp = F)
    legend("topright", c("99%","99.9%"), lty=c(1,1), lwd=c(1,1),col=c("blue", "red"), cex = 0.6, bty = "n")
    dev.off()
    message("-log10(P-values) plotted for ", pairname)
  }
}

message("Finished all tasks!")

message("ENDED> XP-EHH Plot")

#<END>
