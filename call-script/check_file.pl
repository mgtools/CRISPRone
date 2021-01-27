#!/usr/local/perl

#0: fasta 
#1: output fasta
#2: summary 

$seq_check_out="Done";
open (CONT,"$ARGV[0]");
open (NEW,">$ARGV[1]");
%contig_len=();
%contig_des=();
$flag=0;
$n=1;
$short = "NA";
$full = "NA";
$add=0;
while (<CONT>){
    #chomp $_;
    $_ =~ s/\x0d{0,1}\x0a\Z//s;
    if ($_=~/^\>(.+)/){
	$flag++;
	@head=split(/\s+/,$1);
	$short = shift @head;
 	$full = join(" ", @head);
	if (length($short)==0){
	    $seq_check_out="error1: format error";
	    print "$seq_check_out\n";
	    exit;
	}
	if (length($short)>40){
	    $short=substr $short,0,40;
	}
	#$short=$short."_".$n;
	$short=$short; #to check YY June 12, 2018
	$n++;
	$contig_des{$short}=$full;
	$contig_len{$short}=0;
   	print NEW ">$short\n";
    }
    else{
	$_=uc($_);
	$thisline = "";
	if ($_!~/^[A|T|G|C|N]+$/){
	    @letter=split(//,$_);
	    for ($i=0;$i<=@letter-1;$i++){
		 if ($letter[$i]=~/[A|T|G|C|N]/){
		     $thisline = $thisline.$letter[$i]; 
		 }elsif ($letter[$i]=~/[M|R|W|V|H|D|S|Y|B|K|X]/){
		     $thisline = $thisline."N";  
		 }else{
		     $thisline=$thisline.$letter[$i];
		 }
	    }
	}else{
	    $thisline = $_;
	}
	$contig_len{$short} = $contig_len{$short} + length($thisline);
   	print NEW "$thisline\n";
    }
}
close (CONT);
close(NEW);

print "processed $add\n";

if ($flag==0){
    $seq_check_out="error2: not fasta format";
    print "$seq_check_out\n";
    exit;
}else{
    for $key (keys %contig_seq){
	if(($contig_len{$key})==0){
	    $seq_check_out="error3: format error";
	    print "$seq_check_out\n";
	    exit;
	}
    }
}

print "$seq_check_out\n";

open (SUM,">$ARGV[2]");
for $key (keys %contig_len){
   $len=$contig_len{$key};
   $note=$contig_des{$key};
   print SUM "$key\t$len\t$note\n";
}
close(SUM);
