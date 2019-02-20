#Performs haplotype tests per population per chromosome

for pop in POP1 POP2 POP3 POP4
do
	for nchr in {1..22}
	do
		file="chr$nchr/$pop.chr$nchr"
		if [ ! -f "$file.ihs" ]
		then
			echo "$file"
			perl run_conversion_and_parsing.pl $file $nchr
		fi
	done

done
