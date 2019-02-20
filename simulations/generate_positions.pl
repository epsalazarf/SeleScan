#!/usr/bin/env perl

use warnings;
use strict;
use POSIX;
use Cwd;

###

my $in=$ARGV[0];

###

my @positions;

open(IN,$in);
foreach my $line (<IN>){
	chomp($line);
	if($line =~ /positions: (\d+)/){
		@positions=split(" ",$line);
	};
};
close(IN);

#####

open(OUT,">".$in.".pos");
foreach my $pos (@positions){
	print OUT $pos."\n";
};
close(OUT);



