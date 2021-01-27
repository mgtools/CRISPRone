#!/usr/local/perl

$gff_file=$ARGV[1];
open (FILE,"$ARGV[0]");
%begin_pos=();
while (<FILE>){
    if ($_=~/\>(.+)/){
	$head=$1;
	@tmp=split(/\_/,$head);
	$sample=$tmp[0];
	if (@tmp>=5){
	    for ($i=1;$i<=@tmp-4;$i++){
		$sample=$sample."_".$tmp[$i];
	    }
	}
	$key=$sample;
	if (exists $begin_pos{$key}){
	    $begin_pos{$key}=$begin_pos{$key}."\n".$head;
	}else{
	    $begin_pos{$key}=$head;
	}
    }
}
close(FILE);

open (NEW,">$gff_file");
$cds_n=0;
print NEW "##gff-version 3\n";
print NEW "#!gff-spec-version 1.20\n";
print NEW "#!processor Quan Zhang\n";
print NEW "##sequence-region $sample 1 \n";

for $key (keys %begin_pos){ 
    %pos_tmp=();
    @lines=split(/\n/,$begin_pos{$key});
    for ($i=0;$i<=@lines-1;$i++){
	@tmp=split(/\_/,$lines[$i]);
	$begin_pos_pos=@tmp-3;
	$pos_tmp{$tmp[$begin_pos_pos]}=$lines[$i];
    }
    for $key_pos ( sort {$a<=>$b} keys %pos_tmp) { 
	@tmp=split(/\_/,$pos_tmp{$key_pos});
	$end_pos=@tmp-2;
	$dir_pos=@tmp-1;
	print NEW "$key\tRefSeq\tCDS\t$key_pos\t$tmp[$end_pos]\t.\t$tmp[$dir_pos]\t.\tID=$pos_tmp{$key_pos};Name=$pos_tmp{$key_pos};Parent=gene$cds_n\n";
	$cds_n++;
    }
}
print NEW "###\n";
close(NEW);





