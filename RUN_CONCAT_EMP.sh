#Concatenates all empirical data results to single files

#Parameters
pops=("POP1" "POP2" "POP3" "POP4")
pairs=("POP1xPOP2" "POP1xPOP3" "POP1xPOP4" "POP2xPOP3" "POP2xPOP4" "POP3xPOP4")
all="ALL"

echo "STARTED> Concatenate empirical data:"
cd PopGenome/
mkdir results
echo "Concatenating Neutrality..."
bash concat_neut_per_pop.sh "${pops[@]}"
echo "Concatenating all SFS..."
bash concat_sfs_all.sh $all
echo "Concatenating SFS per pop..."
bash concat_sfs_per_pop.sh "${pops[@]}"
echo "Concatenating all Fst..."
bash concat_fst_all.sh $all
echo "Concatenating paired Fst..."
bash concat_fst_per_pair.sh "${pairs[@]}"
cd ../rEHH/
echo "Concatenating iHS..."
bash concat_ihs_per_pop.sh "${pops[@]}"
cd ..
echo "FINISHED> All concatenations for empirical data"
