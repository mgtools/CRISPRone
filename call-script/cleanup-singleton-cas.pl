#!/usr/local/perl

open (FILE,"$ARGV[0]");
@line=<FILE>;
close(FIILE);
@loci_tag=();
@loci_cas=();
$count=0;
for ($i=0;$i<=@line-1;$i++){
    chomp $line[$i];
    @tmp=split(/\s+/,$line[$i]);
    @cas=split(/\:/,$tmp[-2]);
    
    if ($tmp[4]=~/repeat:/){
	$tag="[repeat]";
    }elsif ($tmp[4]=~/antiRepeat/){
	$tag="[antiRepeat]";
    }else{
	$tag=$tmp[4]."[CAS]";
    }
    if ($i==0){
	$loci_cas[$count]=$line[$i];
	$loci_tag[$count]=$tag;
	$pre_end=$tmp[2];
    }else{
	if (($tmp[1]-$pre_end)>=10000){
	    $count++;
	    $loci_tag[$count]=$tag;
	    $loci_cas[$count]=$line[$i];
	    $pre_end=$tmp[2];
	}else{
	    $loci_cas[$count]=$loci_cas[$count]."\n".$line[$i];
	    $loci_tag[$count]=$loci_tag[$count]."-".$tag;
	    $pre_end=$tmp[2];
	}
    }
}


$new=$ARGV[0].".new";
open (NEW,">$new");
for ($i=0;$i<=@loci_tag-1;$i++){
    my $cas_num = () = $loci_tag[$i] =~ /\[CAS\]/g;  # $c is now 3
    #$total_num=@cas_num;
    @tmp=split(/\n/,$loci_cas[$i]);
    if (@tmp==1 && $cas_num==1){
    }else{
	print NEW "$loci_cas[$i]\n";
    }
#    if ($cas_num!=1){
#	print NEW "$loci_cas[$i]\n";
#    }else{
#	@tmp=split(/\n/,$loci_cas[$i]);
#	for ($j=0;$j<=@tmp-1;$j++){
#	    if ($tmp[$j]=~/repeat/ || $tmp[$j]=~/antiRepeat/){
#		print NEW "$tmp[$j]\n";
#	    }
#	}
#    }
}
close(NEW);

`mv $ARGV[0] $ARGV[0].old`;
`mv $new $ARGV[0]`;

