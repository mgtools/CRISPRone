#!/usr/local/perl                                                                                                                                                                    
use List::Util qw(min max);
use File::Basename;
my $dirs = dirname($0); #July 2017, YY 
my $crisprone="$dirs/../";
#0: repeat spacer file
#1: summery filter
#2: output
open (CONT,"$ARGV[0]"); #repeat spacer                                                                                                                                               
%contig=();
while (<CONT>){
    chomp $_;
    if ($_=~/\>(.+)/){
        $head=$1;
	@tmp=split(/\s+/,$head);
	if ($tmp[1]=~/repeat_pos:(\d+)-(\d+)/){
	    $begin=$1;
	    $end=$2;
	}
	if(exists $contig{$tmp[0]}){
	    $contig{$tmp[0]}=$contig{$tmp[0]}.";".$begin."\t".$end."\t.\trepeat\t".$tmp[2];
	}else{
	    $contig{$tmp[0]}=$begin."\t".$end."\t.\trepeat\t".$tmp[2];
	}
        
    }
    else{
        $contig{$tmp[0]}=$contig{$tmp[0]}."\t".$_;
    }
}
close (CONT);

%subtype_name=();
open (FILE,"$crisprone/local/cas-db/cas-CRISPR.GA");
while(<FILE>){
    if ($_=~/^\#/){
	next;
    }
    @tmp=split(/\s+/,$_);

    if ($tmp[2]=~/(.+)\:/){
	$fun_tmp=$1;
	$key=$tmp[0].":".$fun_tmp;
	if ($fun_tmp=~/^cas1$/ || $fun_tmp=~/^cas2$/){
	    $subtype_name{$key}="universal";
	    next;
	}
    }
   
    @each_type=split(/\,/,$tmp[3]);
    my %gen_type=();
    my %spe_type=();
    for ($i=0;$i<=@each_type-1;$i++){
	@tmp_type=split(/\-/,$each_type[$i]);
	if (!exists $gen_type{$tmp_type[1]}){
	    $gen_type{$tmp_type[1]}=1;
	}
	if (@tmp_type==3){
	    if (!exists $spe_type{$each_type[$i]}){
		$spe_type{$each_type[$i]}=1;
	    }
	}
	
    }
    if ((keys %gen_type)==1){
	if ((keys %spe_type)==1){
	    $spe_key=(keys %spe_type)[0];
	    $spe_key=~s/CAS\-/Subtype-/g;
	    $subtype_name{$key}=$spe_key;
	}else{
	    $subtype_name{$key}="Type-".(keys %gen_type)[0];
	}
    }else{
	$subtype_name{$key}="other";
    }
}
close(FILE);

open (ORDER,"$ARGV[1]");  #summary-filter
while (<ORDER>){
    chomp $_;
    if ($_!~/^\#/){
	@tmps=split(/\s+/,$_);
	#$orignal_tmps6=$tmps[6];
	
	if ($tmps[6]=~/^unk$/){
	    $type="unknown";
	    
	}else{
	    $type=$subtype_name{$tmps[6]};
	}

	if (exists $contig{$tmps[2]}){
	    $contig{$tmps[2]}=$contig{$tmps[2]}.";".$tmps[3]." ".$tmps[4]." ".$tmps[5]." ".$tmps[6]." ".$type;
	}else{
	    $contig{$tmps[2]}=$tmps[3]." ".$tmps[4]." ".$tmps[5]." ".$tmps[6]." ".$type;
	}

    }
    
}
close(ORDER);

open (NEW,">$ARGV[2]"); #output

for $key (keys %contig){
    print NEW "##$key\n";
    @tmp=split(/\;/,$contig{$key});
    
    my %pos=();

    for ($i=0;$i<=@tmp-1;$i++){
	@infor=split(/\s+/,$tmp[$i]);
	if ($infor[3]=~/^repeat$/){
	    @array=split(/\-/,$infor[4]); 
	    if (exists $pos{$infor[0]}){
		@pre=split(/\s+/,$pos{$infor[0]});
		if ($pre[2]<$infor[1]){
		    $pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3].":".$infor[5]."\t".$infor[4];
		}
	    }else{
		$pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3].":".$infor[5]."\t".$infor[4];
	    }
	}else{
	    $pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3]."\t".$infor[4];
	}
    }
    for $key ( sort {$a<=>$b} keys %pos) {
        print NEW "$pos{$key}\n";
    }


}
close(NEW);
