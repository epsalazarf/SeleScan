# Split VCF into chromosomes and compress them

for nchr in {1..22}

do

  mkdir chr$nchr

  vcftools --gzvcf ../ALL.SS.adna.phased.vcf.gz --keep ../ALL.ids --chr $nchr --mac 2 --max-missing 1.0 --min-alleles 2 --max-alleles 2 --recode --out chr$nchr/ALL.chr$nchr

  bgzip chr$nchr/ALL.chr$nchr.recode.vcf
  tabix -p vcf chr$nchr/ALL.chr$nchr.recode.vcf.gz

done
