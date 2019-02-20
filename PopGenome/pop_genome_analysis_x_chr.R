#Activated by RUN_POPGENOME.sh

# <START>

# <DATA>

args <- commandArgs(trailingOnly = TRUE)
nchr <- args[1]

chr <- paste('chr', nchr, sep= "")

lchr <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432)

chr_len <- lchr[as.integer(nchr)]

library("PopGenome", quietly=T)

message("Chromosome ", nchr,": Reading data...")

data <- readVCF(paste(chr,"/ALL.",chr,".recode.vcf.gz",sep=""), 10000, nchr, 1, chr_len)
ids <- read.table("../input/pops.ids", stringsAsFactors = F)
pops <- unique(ids$V1)

# </DATA>

# <OUTGROUP>
#outgroup <- c()
#data <- set.outgroup(data,OG)

# <POPULATIONS>

message("Chromosome ", nchr,": Setting populations...")

#Declare populations to be analyzed
for(p in pops){assign(paste0(p),ids[ids$V1==p,2])}
message("Found populations:")
for (i in pops){message(i,": ",length(get(i)), " samples")}

data <- set.populations(data,mget(pops), diploid=T)

# <SLIDING.WINDOWS>
message("Chromosome ", nchr,": Creating windows")

data.slide <- sliding.window.transform(data,width=25000, jump=12500,type=2,whole.data=TRUE)

# <NEUTRALITY.TESTS>

message("Chromosome ", nchr,": Neutrality tests")

data.slide <- neutrality.stats(data.slide,FAST=TRUE)
neutral_tests <- get.neutrality(data.slide)

message("Chromosome ", nchr,": Generating output")

for (i in 1:length(pops)){
  write.table(neutral_tests[i],file=paste0(chr,"/neutrality.results.",chr,".",pops[i],".txt"))
}

message("Chromosome ", nchr,": Finished neutrality tests")

# MAF & FST for ALL
message("Chromosome ", nchr,": Starting Fst and SFS tests")
# <MAF> (SNP-based)
data <- detail.stats(data, site.FST=TRUE)
write.table(t(data@region.stats@minor.allele.freqs[[1]]),file=paste0(chr,"/SFS.results.",chr,".ALL.txt"))

# <FST> (SNP-based)
write.table(data@region.stats@site.FST[[1]],file=paste0(chr,"/fst.results.",chr,".ALL.txt"))

# MAF & FST per POP
for (i in 1:length(pops)){
  write.table(neutral_tests[i],file=paste0(chr,"/neutrality.results.",chr,".",pops[i],".txt"))
  data <- detail.stats(data,new.populations=list(get(pops[i])))
  write.table(t(data@region.stats@minor.allele.freqs[[1]]),file=paste0(chr,"/SFS.results.",chr,".",pops[i],".txt"))
}

# Paired FST
pairs <- combn(pops, 2)
for (i in 1:length(pairs)){
  data<-detail.stats(data, new.populations=mget(pairs[,i]), site.FST=TRUE)
  write.table(data@region.stats@site.FST[[1]],file=paste0(chr,"/fst.results.",chr,".",pairs[1,i],"x",pairs[2,i],".txt"))
}

message("Chromosome ", nchr,": Finished Fst and SFS tests")

message("Chromosome ", nchr,": Finished all tests")

# <END>
