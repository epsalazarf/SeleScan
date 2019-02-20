#Activated by RUN_PIPELINE_3.sh

cd rEHH

for ((i = $1; i <= $2; i++)); do

	perl run_conversion_and_parsing_2.pl sim_$i 

done

cd ..
