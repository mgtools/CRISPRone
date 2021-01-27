#!/usr/local/perl

#0: sequence file
#1: gff-file 
#2: output-prefix-only


open (FILE,"$ARGV[0]");
%seq=();
%head_detail=();
while (<FILE>){
    if ($_=~/^\>(.+)/){
	$head=$1;
	chomp $head;
	$seq{$head}="";

    }
    else{
	chomp $_;
	$seq{$head}=$seq{$head}.$_;

    }

}
close(FILE);

open (FILE,"$ARGV[1]");
$protein_output=$ARGV[2].".faa";
$gene_output=$ARGV[2].".fna";
$gnum=0;
while (<FILE>){
    chomp $_;
    @tmp=split(/\s+/,$_);
    #if ($tmp[2]!=~/repeat/ || $tmp[2]=~/region/){
    if (!($tmp[2]=~/CDS/)) {
	next;
    }
    if (exists $seq{$tmp[0]}){
	if ($gnum==0) {
		open (NEWP,">$protein_output");
		open (NEWG,">$gene_output");
	}
	$gnum = $gnum + 1;
	$gene=substr $seq{$tmp[0]},($tmp[3]-1),($tmp[4]-$tmp[3]+1);
	#if ($tmp[3]=~/\-/){
	#YY July 27 2018
	if ($tmp[6]=~/\-/){
	    $gene=&comp($gene);
	}
	$protein = "";
	for ($j = 0; $j < (length($gene) - 2); $j = $j + 3){
	    $codon=substr($gene, $j, 3);
	    $aa = &codon2aa($codon);
	    if (($aa eq "*") && ($j < (length($gene) - 3))){
		$aa="U";
	    }
	    if ($aa=~/\?/){
		print "error: noncodon encountered!($cds_id:$tmp[0] $tmp[1] $tmp[2])\n";
		exit;
	    }
	    $protein = $protein.$aa;
	}
	$protein=~s/\*//g;
    	@tmp2=split(/;/,$tmp[8]);
    	@tmp3=split(/=/,$tmp2[3]);
	$name=$tmp3[1];
	if ($name eq "") { $newid = "contig_ID:$tmp[0] Gene_location:$tmp[3]-$tmp[4] direction:$tmp[6] Annotation:$tmp[8]"; }
	else { $newid="$name contig_ID:$tmp[0] Gene_location:$tmp[3]-$tmp[4] direction:$tmp[6] Annotation:$tmp[8]"; }
	print NEWP ">$newid\n";
	print NEWG ">$newid\n";
	my @protein_70 = unpack("(A70)*", $protein);
	my @gene_70 = unpack("(A70)*", $gene);
	for ($line=0;$line<=@protein_70-1;$line++){
	    print NEWP "$protein_70[$line]\n";
	}
	for ($line=0;$line<=@gene_70-1;$line++){
	    print NEWG "$gene_70[$line]\n";
	}
    }else{
	print "Cannot find $tmp[0]\n";
    }
}
close(FILE);
if($gnum>0) {
	close(NEWG);
	close(NEWP);
}

sub codon2aa{
    my($codon)=@_;
    $codon=uc $codon;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G','GCN'=>'A','CGN'=>'R','MGR'=>'R','AAY'=>'N','GAY'=>'D','TGY'=>'C','CAR'=>'Q','GAR'=>'E','GGN'=>'G','CAY'=>'H','ATH'=>'I','YTR'=>'L','CTN'=>'L','AAR'=>'K','TTY'=>'F','CCN'=>'P','TCN'=>'S','AGY'=>'S','ACN'=>'T','TAY'=>'Y','GTN'=>'V','TAR'=>'*','TRA'=>'*');
    if(exists $g{$codon})
    {
        return $g{$codon};
    }
    else
    {
        return "?";
    }
}

sub comp{
    my ($str_nu)=@_;
    #print "$str_nu\n";                                                                                                                                                          \
                                                                                                                                                                                  
    @nuclo=split(//,$str_nu);
    for ($n=0;$n<=@nuclo-1;$n++){
        if ($nuclo[$n]=~/T/){
            $nuclo[$n]=~ s/T/A/;
        }
        elsif ($nuclo[$n]=~/A/){
            $nuclo[$n]=~ s/A/T/;
	}
        elsif ($nuclo[$n]=~/C/){
            $nuclo[$n]=~ s/C/G/;
        }
        elsif ($nuclo[$n]=~/G/){
            $nuclo[$n]=~ s/G/C/;
        }
    }
    @re=reverse(@nuclo);
    $re_str=join("",@re);
    return ($re_str);
}
