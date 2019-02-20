#Activated by run_conversion_and_parsing_2.pl


library("rehh", quietly=T)

####

ihh2ihs_neutral<-function(res_ihh,neu_ihh,freqbin=0.05,minmaf=0.05){
  # target ihh
  res_ihh=res_ihh[res_ihh[,3]>minmaf & res_ihh[,3]<(1-minmaf) , ]
  res_ihs=cbind(res_ihh[,1:2],rep(0,nrow(res_ihh)))
  colnames(res_ihs)[3]="iHS" ; rownames(res_ihs)=rownames(res_ihh)
  # neutral ihh
  neu_ihh=neu_ihh[neu_ihh[,3]>minmaf & neu_ihh[,3]<(1-minmaf) , ]
  neu_ihs=cbind(neu_ihh[,1:2],rep(0,nrow(neu_ihh)))
  colnames(neu_ihs)[3]="iHS" ; rownames(neu_ihs)=rownames(neu_ihh)
  if(freqbin>0){
    freq_class=seq(minmaf,1-minmaf,freqbin)
    summary_class=matrix(0,length(freq_class)-1,4) ; colnames(summary_class)=c("Freq Class","Size","Mean iHH","SD iHH")
    ihs_neu=log(neu_ihh[,4]/neu_ihh[,5]) ; ihs_neu[ihs_neu=="Inf" | ihs_neu=="-Inf"]=NA
    ihs_targ=log(res_ihh[,4]/res_ihh[,5]) ; ihs_targ[ihs_targ=="Inf" | ihs_targ=="-Inf"]=NA
     for(c in 1:(length(freq_class)-1)){
       lim_inf=freq_class[c] ; lim_sup=freq_class[c+1]
       mrk_sel=(neu_ihh[,3]>=lim_inf & neu_ihh[,3]<lim_sup)
       mrk_sel_targ=(res_ihh[,3]>=lim_inf & res_ihh[,3]<lim_sup)
       tmp_ihs_neu=ihs_neu[mrk_sel] ; tmp_nmrk_neu=sum(mrk_sel)
       if(tmp_nmrk_neu<10){
          warning(paste("Size of Allele Frequency Class: ",lim_inf,"-",lim_sup," <10 (",tmp_nmrk_neu,"): You should probably increase freqbin\n",sep=""))
       }
       tmp_mean_neu=mean(tmp_ihs_neu,na.rm=T) ; tmp_sd_neu=sd(tmp_ihs_neu,na.rm=T)
       summary_class[c,1]=paste(lim_inf,"-",lim_sup) ; summary_class[c,2]=tmp_nmrk_neu
       summary_class[c,3]=tmp_mean_neu ;  summary_class[c,4]=tmp_sd_neu
       ihs_targ[mrk_sel_targ]=(ihs_targ[mrk_sel_targ]-tmp_mean_neu)/tmp_sd_neu
     }
   }else{
   stop("please, set up freqbin > 0 ")
   }
  res_ihs[,3]=ihs_targ ; tmp_pval=-1*log10(1-2*abs(pnorm(ihs_targ)-0.5))
  tmp_pval2=tmp_pval ; tmp_pval2[tmp_pval2=="Inf"]=NA
  tmp_pval[tmp_pval=="Inf"]=max(tmp_pval2,na.rm=TRUE) + 1
  res_ihs=cbind(res_ihs,tmp_pval) ; colnames(res_ihs)[4]="Pvalue"

  return(list(res.ihs=res_ihs,summary.class=summary_class))
}


####


args <- commandArgs(trailingOnly = TRUE)


load(file=paste("ihh/",args[2],".ihh",sep=""))

freqbin<-unique(sort(ihh_neutral[,3]))

file<-paste(args[1],".",args[2],sep="")

print(file)

#hap<-data2haplohh(hap_file=paste("input/",file,".hap",sep=""),map_file=paste("input/",file,".map",sep=""))

load(file=paste("ihh/",file,".ihh",sep=""))

pdf(paste("plots/",file,".pdf",sep=""))

#res<-scan_hh(hap)

res<-ihh_neutral_temp

ihs<-ihh2ihs_neutral(res,ihh_neutral,freqbin=freqbin[2],minmaf=0.001)
distribplot(ihs$res.ihs[,3])
ihsplot(ihs$res.ihs)
dev.off()

write.table(ihs$res.ihs,file=paste("results/",file,".ihs_temp",sep=""))

print(warnings())

###
