#Concatenates Fst results for each population pair in a single file

for pop in $@
do
	echo $pop
	for nchr in {1..22}
	do
		echo $nchr
		file="chr$nchr/fst.results.chr$nchr.$pop.txt"
		if [ $nchr -eq 1 ]
		then
			echo "CHR POS FST" > results/fst.results.adna.$pop.txt
		fi
		awk -v c=$nchr 'NR>1{print c,$1,$2}' $file >> results/fst.results.adna.$pop.txt
	done
	sed -i 's/"//g' results/fst.results.adna.$pop.txt
done
