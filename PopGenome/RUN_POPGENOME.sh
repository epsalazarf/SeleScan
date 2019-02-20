#Perform neutrality and Fst tests for each chromosome

for i in {1..22}; do
	Rscript pop_genome_analysis_x_chr.R $i
done
