#Performs iHS tests from iHH and plots the results (LONG DURATION)

for ((i = 1; i <= 5000; i=i+100)); do

	n=$((i+99))
	bash RUN_rEHH_2nd_step.sh $i $n

done
