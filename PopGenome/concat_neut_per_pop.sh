#Concatenates neutrality results for each population in a single file

for pop in $@
do
	echo $pop
	for nchr in {1..22}
	do
		echo $nchr
		file="chr$nchr/neutrality.results.chr$nchr.$pop.txt"
		if [ $nchr -eq 1 ]
		then
			echo "CHR start end ss TD FLF FLD" > results/neut.results.adna.$pop.txt
		fi
		awk -v c=$nchr 'NR>1{print c,$1,$3,$6,$5,$8,$9}' $file >> results/neut.results.adna.$pop.txt
	done
	sed -i 's/"//g' results/neut.results.adna.$pop.txt
done
