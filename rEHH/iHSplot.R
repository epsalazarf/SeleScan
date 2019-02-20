# Selection Scan Pipeline: iHS Plot
# By: Pavel Salazar

#Requires: concat_ihs_per_pop.sh output files

#<INPUT>

ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)
suffix <-(".adna.ihs")
outputs <- c("T","P","V") #[T]able, [P]-value plots, [V]enn diagram
plotformat <- ("png") #pdf or png (recommended)

#<START>
message("STARTED> iHS Plot")

# Libraries
library(data.table)

# Functions
antilog10<-function(lx){exp(lx/log(exp(1),10))}

# UCSC Chromosome Palette: (Luminosity 80%)
chrcols <- c("#f4bc7f","#caca88","#cbcb67","#ff8989","#ff7474","#ff77ff","#edbaba","#ffae3e","#f3c000","#cece00","#a6d800","#aef49b","#92d87e","#babaff","#98c4ff","#98cbfe","#00e0e0","#9cd0d0","#faa1ff","#ff92ff","#e3b1ff","#c6c6c6")

#<CORE>
for (p in pops) {
  message("Reading data for ", p,"...")
  file <- (paste(p,suffix, sep= ""))
  results <- data.frame(fread(file, header=T))

  # Remove NAs, absolute iHS and add raw P values
  message("Transforming data for ", p,"...")
  results <- results[complete.cases(results),]
  results$iHS <- abs(results$iHS)
  results$Praw <- antilog10(-results$Pvalue)

  # Top Tables
  if("T" %in% outputs){
    restop1 <- subset(results,results$iHS > quantile(results$iHS,prob=0.99, na.rm = T))
    write.table(restop1,file = paste("results/",p,".ihs.top1.txt",sep = ""), quote = F, sep = "\t", row.names = F)
    message("Wrote Top 1% table for ", p)
    restop01 <- subset(results,results$iHS > quantile(results$iHS,prob=0.999, na.rm = T))
    write.table(restop01,file = paste("results/",p,".ihs.top01.txt",sep = ""), quote = F, sep = "\t", row.names = F)
    message("Wrote Top 0.1% table for ", p)
  }

  # qqman: Manhattan Plot for P-values
  if("P" %in% outputs){
    library(qqman)
    results$Col <- chrcols[results$CHR]
    pvlines <- c(quantile(results$Pvalue,prob=0.99, na.rm = T),quantile(results$Pvalue,prob=0.999, na.rm = T))
    message("Plotting data for ", p,"...")
    message("Format selected: ", plotformat)
    switch(plotformat,
      pdf = pdf(paste("results/",file,".pdf",sep=""), height= 4.5, width= 8),
      png(filename = paste("results/",file,".png",sep=""), width = 1920, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
    )

    # Plot logP values
    manhattan(results, snp = "ID", bp = "POSITION", p= "Pvalue", genomewideline = pvlines[2], suggestiveline = pvlines[1], cex.axis=0.7, cex = 0.7, main= p, col=chrcols, logp = F)
    legend("topright", c("99%","99.9%"), lty=c(1,1), lwd=c(1,1),col=c("blue", "red"), cex = 0.6, bty = "n")
    if (plotformat == "png") {
      dev.off()
      png(paste("results/",file,".qq.png",sep=""), width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
    }
    message("P values plotted for ", p)

    # Q-Q Plot
    qq(results$Praw, main= p, cex.axis=0.7, cex = 0.7)
    dev.off()
    message("Q-Q plotted for ", p)
    message("Finished: ",p)
  }
}

# Venn diagrams
if("V" %in% outputs){
  message("Generating Venn diagram...")
  source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
  for (t in c(1,0.1)){
    tx <- sub(".","",t,fixed = T)
    idlist <- list()
    for (p in pops) {
      temp <- read.table(paste("results/",p,".ihs.top",tx,".txt",sep=""),header=T,colClasses = c(NA,"NULL","NULL","NULL","NULL","NULL"))
      idlist[[p]] <-as.character(temp[,1])
    }
    OLlist <- overLapper(idlist,type="vennsets")
    counts <- list(sapply(OLlist$Venn_List, length),sapply(OLlist$Venn_List, function (x) round(length(x)/dim(OLlist$Intersect_Matrix)[1]*100,2)))

    png(paste("results/ALL.top",tx,".venn.png",sep=""), width = 1080, height = 1080, units = "px", pointsize = 12, bg = "white", res=96, type="cairo")
    vennPlot(counts=counts, mymain= "iHS Shared SNPs", mysub=paste("Top ",t,"%",sep=""), yoffset=c(0.3, -0.2))
    dev.off()
    }
}

message("Finished all tasks!")

message("ENDED> iHS Plot")

#<END>
