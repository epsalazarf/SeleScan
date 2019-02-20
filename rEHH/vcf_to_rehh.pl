#!/usr/bin/perl

use warnings;
use strict;
use POSIX;
use Cwd;

#######

my $file=$ARGV[0];
my $chr=$ARGV[1];

#######

system("vcftools --gzvcf $file.recode.vcf.gz --IMPUTE --out $file.temp");

########

open(IN,"$file.temp.impute.hap");
my @data=<IN>;
my @arr=split(/ /,$data[0]);
close(IN);

my $rows=$#data+1;
my $columns=$#arr+1;

system("transpose -t --fsep \' \' --input ".$rows."x".$columns." $file.temp.impute.hap | sed -e 's/1/2/g' | sed -e 's/0/1/g' > $file.pre.hap");

my $count=1;
open(OUT,">$file.hap");
open(IN,"$file.pre.hap");
foreach my $line (<IN>){
	print OUT "$count $line";
	$count++;
};
close(IN);
close(OUT);

##########

open(OUT,">$file.map");
open(IN,"$file.temp.impute.legend");
my $header=<IN>;
foreach my $line (<IN>){
	chomp($line);
	my @array=split(" ",$line);
	print OUT "$array[0] $chr $array[1] 1 2\n";
};
close(IN);
close(OUT);

###

system("rm $file.pre.hap $file.temp.log $file.temp.impute.legend $file.temp.impute.hap $file.temp.impute.hap.indv");
