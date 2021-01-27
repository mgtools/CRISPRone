#!/usr/bin/perl -Ibin/perl-lib/lib/perl5/ 

#0: input spacer file

use Parallel::ForkManager;
use File::Basename;
my $dirs = dirname($0); #July 2017, YY 
my $crisprone="$dirs/../";

my $max_process = 48;

my $pm = new Parallel::ForkManager($max_process);

#system "date +%T";


$input=$ARGV[0]."/all-spacer.fa";
$db_dir="$crisprone/local/all-mge-db/";
opendir (DB,"$db_dir");
while ($db=readdir(DB)){
    if ($db=~/(.+).fna$/){
	$db_name=$1;
	$pm->start and next; # do the fork                                                                                                                                
	$out_name=$ARGV[0]."/".$db_name.".out";
	`$crisprone/bin/blast+/bin/blastn -db $db_dir/$db -query $input -out $out_name -outfmt 6 -word_size 11 -num_threads 8`;
	$pm->finish;
    }
}

$pm->wait_all_children;


open (FILE,"$crisprone/local/bact-draft-crispr-array.info.rename");
%crispr=();
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if (exists $crispr{$tmp[0]}){
        $crispr{$tmp[0]}=$crispr{$tmp[0]}."\t".$tmp[1]."-".$tmp[2];
    }else{
        $crispr{$tmp[0]}=$tmp[1]."-".$tmp[2];
    }
}
close(FILE);


open (FILE,"$crisprone/local/mge-db/crispr-mge-cdhitest0.8-ann-summary.txt");
%mge_tag=();
while(<FILE>){
    if ($_=~/^\#/){
    }else{
        @tmp=split(/\s+/,$_);
        chomp $_;
        if ($_=~/tag:\s+(.+)/){
            $tag=$1;
            if ($tag=~/(.+)\s+\(/){
                $mge_tag{$tmp[0]}=$1."(HMP)";
            }elsif ($tag=~/UNK/){
                $mge_tag{$tmp[0]}="mge(HMP)";
            }elsif($tag=~/ICE/){
                $mge_tag{$tmp[0]}="ICE(HMP)";
            }
        }
    }
}
close(FILE);

$output=$ARGV[0]."/spacer-mge.result";
open (NEW,">$output");
opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/crispr-mge.out/){
	$sample=$1;
	open(FILE,"$ARGV[0]/$file");
	%gene=();
	#%viral_n=();
	
	while (<FILE>){
	    chomp $_;
	    @tmp=split(/\s+/,$_);
	    @gene_infor=split(/\:/,$tmp[0]);
	    $len_pos=@gene_infor-1;
	  if ($tmp[0]=~/^(.+)_s\d+:\d+$/){
                $id=$1;
                $check=index($tmp[1],$id);
	  }
            if ($check!=-1){
                next;
            }

	    if (($tmp[3]>=($gene_infor[$len_pos]*0.8)) && ($tmp[2]>=0.8)){
		if (exists $mge_tag{$tmp[1]}){
		    print NEW "$_\t$mge_tag{$tmp[1]}\n";
		}else{
		    print NEW "$_\tmge(HMP)\n";
		}
	    }
	}
	close(FILE);
    }
    if ($file=~/all-plasmid.out/){
	$sample=$1;
        open(FILE,"$ARGV[0]/$file");
	while (<FILE>){
	    chomp $_;
            @tmp=split(/\s+/,$_);
	    @gene_infor=split(/\:/,$tmp[0]);
	    $len_pos=@gene_infor-1;
	    if ($tmp[0]=~/^(.+)_s\d+:\d+$/){
                $id=$1;
                $check=index($tmp[1],$id);
            }
            if ($check!=-1){
                next;
            }

	    if (($tmp[3]>=($gene_infor[$len_pos]*0.8)) && ($tmp[2]>=0.8)){
		print NEW "$_\tplasmid(ncbi)\n";
	    }
	}
	close(FILE);
	
    }

    if ($file=~/host-bacteria.out/){
	$sample=$1;
	open(FILE,"$ARGV[0]/$file");
	while (<FILE>){
	    chomp $_;
	    @tmp=split(/\s+/,$_);
	    @gene_infor=split(/\:/,$tmp[0]);
	    $len_pos=@gene_infor-1;
	    if ($tmp[0]=~/^(.+)_s\d+:\d+$/){
                $id=$1;
                $check=index($tmp[1],$id);
            }
            if ($check!=-1){
                next;
            }

	    if (($tmp[3]>=($gene_infor[$len_pos]*0.8)) && ($tmp[2]>=0.8)){
		print NEW "$_\tvirus(host:bacteria)\n";
	    }
	}
	close(FILE);
    }
    

    if ($file=~/host-archea.out/){
	$sample=$1;
	open(FILE,"$ARGV[0]/$file");
	while (<FILE>){
	    chomp $_;
	    @tmp=split(/\s+/,$_);
	    @gene_infor=split(/\:/,$tmp[0]);
	    $len_pos=@gene_infor-1;
	    if ($tmp[0]=~/^(.+)_s\d+:\d+$/){
                $id=$1;
                $check=index($tmp[1],$id);
            }
            if ($check!=-1){
                next;
            }
	 
		if (($tmp[3]>=($gene_infor[$len_pos]*0.8)) && ($tmp[2]>=0.8)){
		    print NEW "$_\tvirus(host:archaea)\n";
		}
	}
	close(FILE);
	
    }
   
    if (($file=~/^bacteria.out/) || ($file=~/^draft.out/) ){
        open (FILE,"$ARGV[0]/$file");
        while (<FILE>){
            @tmp=split(/\s+/,$_);
            $flag=1;
	    chomp $_;
            if (exists $crispr{$tmp[1]}){
                @crispr_pos=split(/\t/,$crispr{$tmp[1]});
		$begin=$tmp[8];
		$end=$tmp[9];
                
                for ($i=0;$i<=@crispr_pos-1;$i++){
                    if ($crispr_pos[$i]=~/(\d+)\-(\d+)/){
			$cbegin=$1;
			$cend=$2;
			if (($begin==$cbegin) || ($end==$cend)|| ($begin>$cbegin && $begin<$cend) || ($end>$cbegin && $end<$cend)){
			    $flag=0;
			    last;
			}
                    }
		}
	    }
	    
	    if ($flag==1){
		@gene_infor=split(/\:/,$tmp[0]);
		$len_pos=@gene_infor-1;
		if (($tmp[3]>=($gene_infor[$len_pos]*0.8)) && ($tmp[2]>=0.8)){
		    if ($file=~/bacteria.out/){
			print NEW "$_\tbacteria-complete\n";
		    }else{
			print NEW "$_\tbacteria-draft\n";
		    }
		}
	    }
	}
        close(FILE);
    }

    
    
}
closedir(DIR);
close(NEW);
`rm $ARGV[0]/*.out`;
#`rm $ARGV[0]/all-spacer.fa`;


