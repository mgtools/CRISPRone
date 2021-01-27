#!/usr/local/perl
use List::Util qw(sum);
#0: plot dir containning .plot.txt
#1: crt file as the name input
use File::Basename;
my $dirs = dirname($0); #July 2017, YY 
my $crisprone="$dirs/../";

$log_file=$ARGV[1];
$log_file=~s/.crt/.infor/g;
%inter=();

open (FILE,"$crisprone/local/cas-db/interference-module.txt"); #add path crisprone, Ye Nov 2016
while (<FILE>){
    @tmp=split(/\s+/,$_);
    $tmp[2]=~s/\://g;
    $key=$tmp[0].":".$tmp[2];
    $inter{$key}=1;
}
close(FILE);

open (NEW,">>$log_file");
opendir (DIR,"$ARGV[0]");
while($file=readdir(DIR)){
    if ($file=~/(.+)\.plot\.txt/){
	$id=$1;
	open (FILE,"$ARGV[0]/$file"); #plot
	@line=<FILE>;
	close(FIILE);
	@loci_pos=();
	@loci_cas=();
	$count=0;
	for ($i=0;$i<=@line-1;$i++){
	    @tmp=split(/\s+/,$line[$i]);
	    @cas=split(/\:/,$tmp[-2]);
	    if ($tmp[4]=~/repeat:/){
		$tag="repeat";
	    }elsif ($tmp[4]=~/antiRepeat/){
		$tag="tracrRNA";
	    }elsif ($tmp[4]=~/^unk$/){
		$tag="UNK:UNK[unk](".$tmp[1]."-".$tmp[2].")";
	    }else{
		$tag=$tmp[4]."[CAS](".$tmp[1]."-".$tmp[2].")";
	    }
	    if ($i==0){
		$loci_pos[$count]=$tmp[1];
		$loci_cas[$count]=$tag;
		$pre_end=$tmp[2];
	    }else{
		if (($tmp[1]-$pre_end)>=10000){
		    $loci_pos[$count]=$loci_pos[$count]."-".$pre_end;
		    $count++;
		    $loci_cas[$count]=$tag;
		    $loci_pos[$count]=$tmp[1];
		    $pre_end=$tmp[2];
		}else{
		    $loci_cas[$count]=$loci_cas[$count]."->".$tag;
		    $pre_end=$tmp[2];
		}
	    }
	}
	$loci_pos[$count]=$loci_pos[$count]."-".$pre_end;

	for ($i=0;$i<=@loci_pos-1;$i++){
	    $flag_cas_inter=0;
	    @check_array_inter=split(/\-\>/,$loci_cas[$i]);
	    for ($e=0;$e<=@check_array_inter-1;$e++){
		if ($check_array_inter[$e]=~/(.+)\[CAS\]/){
		    $family_id=$1;
		    if (exists $inter{$family_id}){
			$flag_cas_inter++;
		    }
		}
	    }
	    my $cas_num = () = $loci_cas[$i] =~ /\[CAS\]/g;  # $c is now 3
	    for ($e=0;$e<=@check_array_inter-1;$e++){
		if ($check_array_inter[$e]=~/(.+)\[.+\]\((.+)\)/){
		    #print NEW "\#$id\:cas_pos:$1\:$2\tsuspicious\n";
		    print NEW "\#$id\:cas_pos:$1\:$2\t$flag_cas_inter\t$cas_num\n";
		}
	    }
	}
    }
}
closedir(DIR);
close(NEW);
