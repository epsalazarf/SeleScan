#!/usr/bin/perl

use warnings;
use strict;
use POSIX;
use Cwd;

#######

my $file=$ARGV[0];
my $chr=$ARGV[1];

system("perl vcf_to_rehh.pl $file $chr");
system("Rscript run_rehh_by_chr.R $file");

open(OUT,">".$file.".ihs");
open(IN,$file.".ihs_temp");
my $header=<IN>;
print OUT "\"ID\" ".$header;
foreach my $line (<IN>){
	chomp($line);
	my @array=split(" ",$line);
	print OUT "$array[0] $chr $array[2] $array[3] $array[4]\n";
};
close(IN);
close(OUT);

system("rm ".$file.".ihs_temp");
