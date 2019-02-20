# Merges all simulations' iHH

cd rEHH

for pop in POP1 POP2 POP3 POP4
do
	Rscript merge_iHH.R $pop 10000
done

cd ..
