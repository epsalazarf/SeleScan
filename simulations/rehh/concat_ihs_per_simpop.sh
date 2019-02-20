#Concatenates iHS results for each population in a single file
mkdir concat
for pop in $@
do
	echo $pop
	for nsim in {1..10000}
	do
		if [ $(( $nsim % 100 )) -eq 0 ] ; then
			echo $nsim
		fi
		file="results/sim_$nsim.$pop.ihs"
		if [ $nsim -eq 1 ]
		then
			echo $nsim
			head -1 $file > concat/$pop.sims.ihs
		fi
			awk 'NR>1{print}' $file >> concat/$pop.sims.ihs
	done
	sed -i 's/"//g' concat/$pop.sims.ihs
done
