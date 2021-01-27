#!/usr/local/perl

open (FILE,"$ARGV[0]");
$new=$ARGV[0].".modify";
open (NEW,">$new");
while (<FILE>){
    if ($_=~/\#/){
	print NEW "$_";
    }else{
	@tmp=split(/\s+/,$_);
	$tmp[1]=~s/\-/$tmp[0]/g;
	for ($i=0;$i<=@tmp-1;$i++){
	    print NEW "$tmp[$i]\t";
	}
	print NEW "\n";
    }
}

close(FILE);
close(NEW);

#`rm $ARGV[0]`;
`mv $new $ARGV[0]`;
