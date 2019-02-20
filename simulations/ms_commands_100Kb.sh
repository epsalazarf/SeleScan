#!/bin/bash

#Activated by RUN_ms.sh

#Parameters:
# Loci= 2n (n diploid samples)
# Mutation Rate (μ) = 1.25e-8
# Population Size (N) = ?
# Sequence Length= 100,000 bp
# Theta (θ = 4 ∗ N ∗ μ ∗ Sequence Length) = ?
# Recombination rate range (r=min:max)= 1e-9, 1e-8
# Rho_min (ρ= 4 ∗ N ∗ r min ∗ Sequence Length)= ?
# Rho_max (ρ= 4 ∗ N ∗ r max ∗ Sequence Length)= ?
# Populations (number)= POP1 (?), POP2 (?), POP3 (?), POP4 (?)

# Command
ms [loci] 1 -t [theta] -r $(( ( RANDOM % [rho_max - rho min + 1])  + rho_min )) [seq length +1] -I 4 [pop1] [pop2] [pop3] [pop4] > output.ms

perl add_ancestral.pl input_data/sim_$1.ms
perl generate_positions.pl input_data/sim_$1.ms
