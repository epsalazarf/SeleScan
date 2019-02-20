#1. Create simulations
#2. Performs neutrality and Fst tests with PopGenome
#3. Transforms data to rEHH input and calculates iHH

for ((i = 1; i <= 10000; i=i+100)); do

	n=$((i+99))
	echo $i $n

	bash RUN_first_step.sh $i $n

done
