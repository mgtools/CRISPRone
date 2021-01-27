#!/usr/local/perl                                                                                                                                                                    
use List::Util qw(min max);
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

=head
open (FILE,"local/cas-db/subtype-cas.update.txt");
%subtype_hmm=();
%subtype_name=();
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if (length($tmp[2])>1){
        $tmp[2]=lc($tmp[2]);
        $subtype_hmm{$tmp[2]}=$tmp[1];
    }else{
        $tmp[0]=lc($tmp[0]);
        $subtype_name{$tmp[0]}=$tmp[1];
    }
}
close(FILE);
=cut

%subtype_name=();
open (FILE,"local/cas-db/cas-CRISPR2015.GA");
while(<FILE>){
    if ($_=~/^\#/){
	next;
    }
    @tmp=split(/\s+/,$_);
    if ($tmp[2]=~/(.+)\:/){
	$fun_tmp=$1;
	$key=$tmp[0].":".$tmp[2];
	if ($fun_tmp=~/^cas1$/ || $fun_tmp=~/^cas2$/){
	    $subtype_name{$key}="universal";
	    next;
	}
    }
   
    @each_type=split(/\,/$tmp[3]);
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
	    $subtype_name{$key}=(keys %spe_type)[0];
	}else{
	    $subtype_name{$key}=(keys %gen_type)[0];
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
	$orignal_tmps6=$tmps[6];
	
	if ($tmps[6]=~/^unk$/){
	    $type="unknown";
	    
	}else{
	    $tmps[6]=lc($tmps[6]);
            @cas_ann=split(/\:/,$tmps[6]);

	    if (exists $subtype_hmm{$cas_ann[0]}){
		$type=$subtype_hmm{$cas_ann[0]};
		
	    }else{
		if ($cas_ann[1]=~/cxxc_cxxc/){
		    $cas="cxxc-cxxc";
		    $type=$subtype_name{$cas};
		}elsif($cas_ann[1]=~/crispr\_(.+)/){
		    $cas=$1;
		    if (exists $subtype_name{$cas}){
			$type=$subtype_name{$cas};
		    }else{
			$type="other";
		    }
		}elsif (($cas_ann[1]=~/^cas_cas1$/) || ($cas_ann[1]=~/^cas_cas2$/) || ($cas_ann[1]=~/^cas1$/) ||($cas_ann[1]=~/^cas2$/)){
		    $type="universal";
		}else{
		    if ($cas_ann[1]=~/_/){
			@all_name=split(/\_/,$cas_ann[1]);
			$name_flag=0;
			for ($nu=0;$nu<=@all_name-1;$nu++){
			    if (exists $subtype_name{$all_name[$nu]}){
				$type=$subtype_name{$all_name[$nu]};
				$name_flag=1;
				last;
			    }
			}
			if ($name_flag==0){
			    $type="other";
			}
		    }else{
			if (exists $subtype_name{$cas_ann[1]}){
			    $type=$subtype_name{$cas_ann[1]};
			}else{
			    $type="other";
			}
		    }
		    
		}
	    }
	  
	}
	if (exists $contig{$tmps[2]}){
	    $contig{$tmps[2]}=$contig{$tmps[2]}.";".$tmps[3]." ".$tmps[4]." ".$tmps[5]." ".$orignal_tmps6." ".$type;
	}else{
	    $contig{$tmps[2]}=$tmps[3]." ".$tmps[4]." ".$tmps[5]." ".$orignal_tmps6." ".$type;
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
	    #if ($array[-1]=~/r(\d+)/){
	#    $last_num=$1;
#		 if ($last_num>=4){
	    if (exists $pos{$infor[0]}){
		@pre=split(/\s+/,$pos{$infor[0]});
		if ($pre[2]<$infor[1]){
		    $pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3].":".$infor[5]."\t".$infor[4];
		}
	    }else{
		$pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3].":".$infor[5]."\t".$infor[4];
	    }
#		 }
	    # }

#	    $pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3].":".$infor[5]."\t".$infor[4];
	    #print NEW "$key\t$infor[0]\t$infor[1]\t$infor[2]\t$infor[3]\:$infor[5]\t$infor[4]\n";
	}else{
	    $pos{$infor[0]}=$key."\t".$infor[0]."\t".$infor[1]."\t".$infor[2]."\t".$infor[3]."\t".$infor[4];
	    #print NEW "$key\t$infor[0]\t$infor[1]\t$infor[2]\t$infor[3]\t$infor[4]\n";
	}
    }
    for $key ( sort {$a<=>$b} keys %pos) {
        print NEW "$pos{$key}\n";
    }


}
close(NEW);
