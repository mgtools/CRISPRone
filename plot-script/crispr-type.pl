#!/usr/local/perl
#$crisprone="/var/www/html/CRISPRone";
#open (FILE,"$crisprone/local/cas-db/subtype.color");
open (FILE, $ARGV[0]);
%cas_col=();
while (<FILE>){
    @tmp=split(/\s+/,$_);
    $cas_col{$tmp[0]}=$tmp[2];
}
close(FILE);

@begin_array=();
@end_array=();
@dir_array=();
@ann_array=();
@col_array=();
$overall_begin="[";
$overall_end="[";

$repeat_n=0;
$rna_n=0;
$unk_n=0;
$typeI=0;
$typeII=0;
$typeIII=0;
$typeIV=0;
$typeV=0;
$typeVI=0;
$all=0;
$repeat_poss=0;

@crisp_array=();
$cri_n=0;

$loci_n=0;
$last_pos=-1000000000000000000;
%crispr_array_locus=();

open (FILE,"$ARGV[1]");
while (<FILE>){
    $all++;
    if ($_=~/antiRepeat/){
	$rna_n++;
    }
    @tmp=split(/\s+/,$_);
    if ($tmp[5]=~/ype-VI/) { $typeVI++; }
    elsif ($tmp[5]=~/ype-V/) { $typeV++; }
    elsif ($tmp[5]=~/ype-IV/) { $typeIV++; }
    elsif ($tmp[5]=~/ype-III/){ $typeIII++; }
    elsif ($tmp[5]=~/ype-II/){ $typeII++; }
    elsif ($tmp[5]=~/ype-I/){ $typeI++; }
    elsif ($_=~/unknown/){ $unk_n++; }

    $overall_begin=$overall_begin.$tmp[1].",";
    $overall_end=$overall_end.$tmp[2].",";
    if ($tmp[4]=~/repeat:(.+)/){
	$repeat_n++;
	$ann="repeat";
	$repeat_seq=$1;
	if ($tmp[4]=~/possible/){
	    $col="possible CRISPR array";
	    $repeat_poss++;
	}
        else{
	    @repeat_tmp=split(/\-/,$tmp[5]);
	    $re_n=0;
	    @spacer_tmp=();
	    for ($re=0;$re<=@repeat_tmp-1;$re++){
		if ($repeat_tmp[$re]=~/s(\d+)\:/){
		    $spacer_tmp[$re_n]=$1;
		    $re_n++;
		}
	    }
	    $last_spacer=@spacer_tmp-1;
	    $crisp_array[$cri_n]=$tmp[1]."-".$tmp[2];	    
	    $col="repeat number ".($spacer_tmp[$last_spacer]-$spacer_tmp[0]+2) ."; spacers ".$crisp_array[$cri_n]; 
	    
	    $crispr_array_locus{$crisp_array[$cri_n]}="";
	    $cri_n++;
	}
    }
    elsif ($tmp[4]=~/antiRepeat/){
            $ann="antiRepeat";
            $col=$tmp[5];
    }
    else{
	if (exists $cas_col{$tmp[5]}){
	    $ann=$tmp[4];
	    $col=$cas_col{$tmp[5]};
	}
        else{ $ann="other"; $col="grey"; }
    }
    @new_ann=split(/\:/,$ann);
    if (length($new_ann[1])==0) { $new_ann[1]=$new_ann[0]; }
    if ($last_pos==-1000000000000000000){
        $begin_array[$loci_n]=$tmp[1];
        $end_array[$loci_n]=$tmp[2];
        $dir_array[$loci_n]="\"".$tmp[3]."\"";
	$col_array[$loci_n]="\"".$col."\"";
	$ann_array[$loci_n]="\"".$new_ann[1]."\"";
    }
    else{
	if (($tmp[1]-$last_pos)>=10000){
	    $loci_n++;
	    $begin_array[$loci_n]=$tmp[1];
	    $end_array[$loci_n]=$tmp[2];
	    $dir_array[$loci_n]="\"".$tmp[3]."\"";
	    $col_array[$loci_n]="\"".$col."\"";
	    $ann_array[$loci_n]="\"".$new_ann[1]."\"";
	}
        else{
	    $begin_array[$loci_n]=$begin_array[$loci_n].",".$tmp[1];
	    $end_array[$loci_n]=$end_array[$loci_n].",".$tmp[2];
	    $dir_array[$loci_n]=$dir_array[$loci_n].",\"".$tmp[3]."\"";
	    $col_array[$loci_n]=$col_array[$loci_n].",\"".$col."\"";
	    $ann_array[$loci_n]=$ann_array[$loci_n].",\"".$new_ann[1]."\"";
	}
    }
    #if ($col=~/;\s+spacers\s+(s\d+\-s\d+)/){
    if ($col=~/;\s+spacers\s+(.+\-.+)/) { 
	$crispr_array_locus{$1}=$loci_n;
    }
    $last_pos=$tmp[2];
}
close(FILE);

$overall_begin =~ s/(.+),$/$1\]/;
$overall_end =~ s/(.+),$/$1\]/;

$cas=$all-$repeat_n-$unk_n-$rna_n;
$type_str="NA";
@types=();
if($typeI>=1){ push(@types, "I"); }
if($typeII>=1){ push(@types, "II"); }
if($typeIII>=1){ push(@types, "III"); }
if($typeIV>=1){ push(@types, "IV"); }
if($typeV>=1){ push(@types, "V"); }
if($typeVI>=1){ push(@types, "VI"); }
if(@types > 0) { $type_str = "(type ".join(", ", @types).")"; }

$print_repeat=$repeat_n-$repeat_poss;
print "$repeat_n\n";
print "$cas\n";
print "$type_str\n";
if (@crisp_array>0){
   for ($i=0;$i<=@crisp_array-1;$i++){ print "$crisp_array[$i]\:$crispr_array_locus{$crisp_array[$i]}\t"; }
}
else{ print "NA"; }
print "\n";

print "$overall_begin\n";
print "$overall_end\n";

$locus_max=0;
for ($i=0;$i<=@begin_array-1;$i++){
    print "$begin_array[$i]\*";
    @locus_begin=split(/\,/,$begin_array[$i]);
    @locus_end=split(/\,/,$end_array[$i]);
    $locus_max_tmp=$locus_end[-1]-$locus_begin[0];
    if ($locus_max_tmp>$locus_max){
	$locus_max=$locus_max_tmp;
    }
}
print "\n";

for ($i=0;$i<=@end_array-1;$i++){ print "$end_array[$i]\*"; }
print "\n";
for ($i=0;$i<=@dir_array-1;$i++){ print "$dir_array[$i]\*"; }
print "\n";
for ($i=0;$i<=@ann_array-1;$i++){ print "$ann_array[$i]\*"; }
print "\n";
for ($i=0;$i<=@col_array-1;$i++){ print "$col_array[$i]\*"; }
print "\n";

print "locus_max:$locus_max\n";
