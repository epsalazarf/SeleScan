#!/usr/bin/perl

use warnings;
use strict;
use POSIX;
use Cwd;

#######

my $sim=$ARGV[0];

my @pops=("POP1","POP2","POP3","POP4");

foreach my $pop (@pops){

	my $file=$sim.".".$pop;

	system("Rscript run_rehh_iHS_x_sim.R $sim $pop");

	open(OUT,">results/".$file.".ihs");
	open(IN,"results/".$file.".ihs_temp");
	my $header=<IN>;
	print OUT "\"ID\" ".$header;
	foreach my $line (<IN>){
		chomp($line);
		my @array=split(" ",$line);
		print OUT "$array[0] 1 $array[2] $array[3] $array[4]\n";
	};
	close(IN);
	close(OUT);

	system("rm results/".$file.".ihs_temp");

};
