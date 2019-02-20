# Split VCF into populations and chromosomes

for pop in POP1 POP2 POP3 POP4
do
  for nchr in {1..22}
  do
    mkdir chr$nchr

    vcftools --gzvcf ../ALL.SS.adna.phased.vcf.gz --keep ../$pop.ids --chr $nchr --mac 2 --max-missing 1.0 --min-alleles 2 --max-alleles 2 --recode --out chr$nchr/$pop.chr$nchr

    bgzip chr$nchr/$pop.chr$nchr.recode.vcf
    tabix -p vcf chr$nchr/$pop.chr$nchr.recode.vcf.gz

  done

done
