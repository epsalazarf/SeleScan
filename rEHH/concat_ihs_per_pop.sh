#Concatenates iHS results for each population in a single file
mkdir results
for pop in $@
do
	echo $pop
	for nchr in {1..22}
	do
		echo $nchr
		file="chr$nchr/$pop.chr$nchr.ihs"
		if [ $nchr -eq 1 ]
		then
			cat $file > results/$pop.adna.ihs
		else
			awk 'NR>1{print}' $file >> results/$pop.adna.ihs
		fi
	done
	sed -i 's/"//g' results/$pop.adna.ihs
done
