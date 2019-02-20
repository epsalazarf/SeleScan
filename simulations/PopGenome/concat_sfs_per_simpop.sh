#Concatenates SFS results for each population in a single file

for pop in $@
do
	echo $pop
	for nsim in {1..10000}
	do
		if [ $(( $nsim % 100 )) -eq 0 ] ; then
			echo $nsim
		fi
		file="results/SFS.sim_$nsim.$pop.txt"
		if [ $nsim -eq 1 ]
		then
			echo $nsim
			echo "SIM POS SFS" > concat/SFS.sims.$pop.txt
		fi
		awk -v s=$nsim 'NR>1{print s,$0}' $file >> concat/SFS.sims.$pop.txt
	done
	sed -i 's/"//g' concat/SFS.sims.$pop.txt
done
