#!/usr/bin/perl

use warnings;
use strict;
use POSIX;
use Cwd;

#######

my $file=$ARGV[0];
my $length=$ARGV[1];

#######


#ALL

my %POP1=("1"=>0,"2"=>0);

my %POP2=("3"=>0,"4"=>0);

my %POP3=("5"=>0,"6"=>0);

my %POP4=("7"=>0,"8"=>0);

########

my %data;
my @positions;
my $read=0;
my $count=1;

my %present_pos;

open(IN,"../input_data/".$file.".ms");
foreach my $line (<IN>){
	chomp($line);
	if($read==1){

		$line =~ s/1/2/g;
		$line =~ s/0/1/g;

		my @haplos=split("",$line);
		$data{$count}=\@haplos;
		$count++;

		if($count==115){
			$read=0;
		};
	};
	if($line =~ /positions: (.+)/){
		my @positions_array=split(" ",$1);

		foreach my $val (@positions_array){

			my $pos=int($val*$length);

			if(exists($present_pos{$pos})){
				$pos++;
			}elsif(exists($present_pos{$pos})){
				$pos++;
			}elsif(exists($present_pos{$pos})){
				$pos++;
			}elsif(exists($present_pos{$pos})){
				$pos++;
			}elsif(exists($present_pos{$pos})){
				$pos++;
			}elsif(exists($present_pos{$pos})){
				$pos++;
			};
			$present_pos{$pos}=0;

			push(@positions,$pos);

		};

		$read=1;
	};
};
close(IN);

####### POP1

open(OUT,">input/".$file.".POP1.map");
foreach my $val (@positions){
	print OUT "chr1_".$val." 1 ".$val." 1 2\n";
};
close(OUT);

open(OUT,">input/".$file.".POP1.hap");
foreach my $id (sort{$a<=>$b}(keys(%data))){
	if(exists($POP1{$id})){
		print OUT $id." @{$data{$id}}\n";
	};
};
close(OUT);

####### POP2

open(OUT,">input/".$file.".POP2.map");
foreach my $val (@positions){
	print OUT "chr1_".$val." 1 ".$val." 1 2\n";
};
close(OUT);

open(OUT,">input/".$file.".POP2.hap");
foreach my $id (sort{$a<=>$b}(keys(%data))){
	if(exists($POP2{$id})){
		print OUT $id." @{$data{$id}}\n";
	};
};
close(OUT);

####### POP3

open(OUT,">input/".$file.".POP3.map");
foreach my $val (@positions){
	print OUT "chr1_".$val." 1 ".$val." 1 2\n";
};
close(OUT);

open(OUT,">input/".$file.".POP3.hap");
foreach my $id (sort{$a<=>$b}(keys(%data))){
	if(exists($POP3{$id})){
		print OUT $id." @{$data{$id}}\n";
	};
};
close(OUT);

####### POP4

open(OUT,">input/".$file.".POP4.map");
foreach my $val (@positions){
	print OUT "chr1_".$val." 1 ".$val." 1 2\n";
};
close(OUT);

open(OUT,">input/".$file.".POP4.hap");
foreach my $id (sort{$a<=>$b}(keys(%data))){
	if(exists($POP4{$id})){
		print OUT $id." @{$data{$id}}\n";
	};
};
close(OUT);
