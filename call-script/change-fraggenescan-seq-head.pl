#!/usr/local/perl

open (FILE,"$ARGV[0]");
open (NEW,">>$ARGV[1]");
while (<FILE>){
    if ($_=~/>(.+)/){
	$head=$1;
	$head=~s/\s+//;
	print NEW ">$head\n";
    }else{
	print NEW "$_";
    }
}
close(FILE);
`mv $ARGV[1] $ARGV[0]`;
