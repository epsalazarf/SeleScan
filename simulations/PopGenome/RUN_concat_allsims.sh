#Concatenates all PopGenome results to directory

mkdir concat
echo "Concatenating Neutrality..."
bash concat_neut_per_simpop.sh
echo "Concatenating all SFS..."
bash concat_sfs_all_sim.sh
echo "Concatenating SFS per pop..."
bash concat_sfs_per_simpop.sh
echo "Concatenating all Fst..."
bash concat_all_simfst.sh
echo "Concatenating paired Fst..."
bash concat_fst_per_simpair.sh
echo "All concatenations finished"
