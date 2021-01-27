#!/usr/local/perl

#0: crt file
#1: output file

open (FILE,"$ARGV[0]"); #.crt file
%hash_repeat=();
while (<FILE>){
    if ($_=~/^SEQ:/){
	@tmp=split(/\s+/,$_);
	#@id=split(/\|/,$tmp[1]);
	#$sample=$id[3];
	$sample=$tmp[1];
    }
    @col=split(/\s+/,$_);
    if (@col==2){
	if ($col[0]=~/^\d+/){
	    $n++;
	    $num_repeat="r".$n;
	    $key=$key."-".$num_repeat;
	    
	}
    }
    if ($_=~/CRISPR\s+\d+\s+Range:\s+(.+)\s+\-\s+(.+)/){
	$start=$1;
	$end=$2;
	if ($start<0){
	    $start=0;
	}
	$key=$sample." repeat_pos:".$start."-".$end." ";
	$n=0;
	next;
    }
    if ($_=~/^\-?\d+\s+.+\s+(.+)\s+\[\s+\d+\,\s+\d+\s+\]/){
	$spacer=$1;
	$n++;
	$num_repeat="r".$n;
	if ($n==1){
	    $key=$key.$num_repeat."-".$spacer;
	}else{
	    $key=$key."-".$num_repeat."-".$spacer;
	}
	
    }
    
    
    if ($_=~/Repeat-con\s+(.+)/){
	$repeat=$1;
#	$hash_repeat{$key}=$repeat;
    }
    if ($_=~/^Repeats:\s+(\d+)\s+Average\s+Length:\s+\d+/){
	#if ($1>=4){
	$hash_repeat{$key}=$repeat; 
	#}
    }
}

open (NEW,">$ARGV[1]");
for $key (keys %hash_repeat){
    print NEW ">$key\n";
    print NEW "$hash_repeat{$key}\n";
}
close(NEW);


