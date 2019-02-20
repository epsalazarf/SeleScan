#Activated by RUN_first_step.sh

for ((i = $1; i <= $2; i++)); do

	bash ms_commands_100Kb.sh $i

done
