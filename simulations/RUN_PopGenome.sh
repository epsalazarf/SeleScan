#Activated by RUN_PIPELINE_1.sh

cd PopGenome

for ((i = $1; i <= $2; i++)); do

	if [ ! -f "results/neutrality.sim_$i.txt" ];
	then
		Rscript pop_genome_analysis_simulations.R sim_$i 100000
	fi

done

cd ..
