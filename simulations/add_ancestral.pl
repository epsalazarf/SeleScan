#!/usr/bin/env perl

use warnings;
use strict;
use POSIX;
use Cwd;

###

my $in=$ARGV[0];

###

my $num_sites;

open(IN,$in);
foreach my $line (<IN>){
	chomp($line);

	if($line =~ /segsites: (\d+)/){
		$num_sites=$1;
	};
};
close(IN);

#####

open(OUT,">>".$in);
for(my $i=1;$i<=$num_sites;$i++){
	print OUT "0";
};
print OUT "\n";
close(OUT);

system("sed -i \'1s\/.*\/ms 9 1\/\' $ARGV[0]");
