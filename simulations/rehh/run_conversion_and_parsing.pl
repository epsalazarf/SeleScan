#!/usr/bin/perl

use warnings;
use strict;
use POSIX;
use Cwd;

#######

my $sim=$ARGV[0];
my $length=$ARGV[1];

my @pops=("POP1","POP2","POP3","POP4");

system("perl ms_to_rehh.pl $sim $length");

foreach my $pop (@pops){

	my $file=$sim.".".$pop;

	system("Rscript run_rehh_iHH_x_sim.R $file");

};
