#Concatenates neutrality results for each population in a single file

for pop in $@
do
	echo $pop
	for nsim in {1..10000}
	do
		if [ $(( $nsim % 100 )) -eq 0 ] ; then
			echo $nsim
		fi
		file="results/neutrality.sim_$nsim.$pop.txt"
		if [ $nsim -eq 1 ]
		then
			echo $nsim
			echo "SIM start end ss TD FLF FLD FWH ZE" > concat/neut.sims.$pop.txt
		fi
		awk -v s=$nsim 'NR>1{print s,$1,$3,$6,$5,$8,$9,$11,$12}' $file >> concat/neut.sims.$pop.txt
	done
	sed -i 's/"//g' concat/neut.sims.$pop.txt
done
