#Concatenates SFS results for each population in a single file

for pop in $@
do
	echo $pop
	for nchr in {1..22}
	do
		echo $nchr
		file="chr$nchr/SFS.results.chr$nchr.$pop.txt"
		if [ $nchr -eq 1 ]
		then
			echo "CHR POS SFS" > results/SFS.results.adna.$pop.txt
		fi
		cut -f2 -d " " $file | sed '1d' > tmpa
		cut -f1 -d " " chr$nchr/fst.results.chr$nchr.ALL.txt | sed '1d' > tmpb
		paste -d, tmpb tmpa > tmpc
		awk -v c=$nchr '{print c,$0}' tmpc >> results/SFS.results.adna.$pop.txt
	done
	sed -i 's/"//g' results/SFS.results.adna.$pop.txt
	sed -i 's/,/ /g' results/SFS.results.adna.$pop.txt
	rm tmp*
done
