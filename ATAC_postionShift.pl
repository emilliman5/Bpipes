#!/usr/local/bin/perl

use strict;

open (SAMIN, "<$ARGV[0]") || die "Could not open $ARGV[0]\n";
open (SAMOUT,">>$ARGV[1]");

while(<SAMIN>){
	chomp;
	
	if ($_=~/^@/){
		print SAMOUT $_."\n";
		}
		
	else{
		
		my @line=split($_, "\t");
		
		if($line[1]==163 | $line[1]==147 ){
			$line[3]=$line[3]-4;
			$line[7]=$line[7]+5;
			$line[8]=$line[8]+1;
		}
		
		elsif($line[1]==83 | $line[1]==99){
			$line[3]=$line[3]+5;
			$line[7]=$line[7]-4;
			$line[8]=$line[8]+1;
			}
		print SAMOUT join("\t",@line)."\n";
		}
	}