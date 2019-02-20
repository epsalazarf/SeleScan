#Activated by RUN_first_step.sh

cd rEHH

for ((i = $1; i <= $2; i++)); do

	if [ ! -f "ihh/sim_$i.ihh" ];
	then
		perl run_conversion_and_parsing.pl sim_$i 100000
	fi


done

cd ..
