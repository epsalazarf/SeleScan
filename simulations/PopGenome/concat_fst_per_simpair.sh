#Concatenates Fst results for each population pair in a single file

for pop in $@
do
	echo $pop
	for nsim in {1..10000}
	do
		if [ $(( $nsim % 100 )) -eq 0 ] ; then
			echo $nsim
		fi
		file="results/fst.results.sim_$nsim.$pop.txt"
		if [ $nsim -eq 1 ]
		then
			echo $nsim
			echo "SIM POS FST" > concat/fst.sims.$pop.txt
		fi
		awk -v s=$nsim 'NR>1{print s,$0}' $file >> concat/fst.sims.$pop.txt
	done
	sed -i 's/"//g' concat/fst.sims.$pop.txt
done
