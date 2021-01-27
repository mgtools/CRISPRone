#!/usr/bin/perl -Ibin/perl-lib/lib/perl5/                                                


#spacer folder 
#repeat-spacer file

use Parallel::ForkManager;
use File::Basename;
my $dirs = dirname($0); #July 2017, YY 
my $crisprone="$dirs/../";

my $max_process = 48;

my $pm = new Parallel::ForkManager($max_process);
%repeat=();
open (FILE,"$ARGV[1]");
$new_repeat=$ARGV[0]."/all-repeat.fa";
$repeat_n=0;
while (<FILE>){
    if ($_=~/>(.+)/){
        $repeat_n++;
    }else{
        chomp $_;
        $repeat{$repeat_n}=$_;
    }
}

open (NEW,">$new_repeat");
for $key (keys %repeat){
    $len=length($repeat{$key});
    print NEW ">$key:$len\n";
    print NEW "$repeat{$key}\n";
}
close(NEW);

$pm->wait_all_children;

$db_dir="$crisprone/local/all-mge-db/";
opendir (DB,"$db_dir");
while ($db=readdir(DB)){
    if ($db=~/(.+).fna$/){
        $db_name=$1;
        $pm->start and next; # do the fork                                          
        $out_name=$ARGV[0]."/".$db_name.".repeat.out";
        `$crisprone/bin/blast+/bin/blastn -db $db_dir/$db -query $new_repeat -out $out_name -outfmt 6 -word_size 11 -num_threads 8`;
        $pm->finish;
    }
}
$pm->wait_all_children;


opendir (DIR,"$ARGV[0]");
%crispr_hit=();
while ($repeat=readdir(DIR)){
    if ($repeat=~/.repeat.out/){
        open (FILE,"$ARGV[0]/$repeat");
        while (<FILE>){
            @tmp=split(/\s+/,$_);
            @len=split(/\:/,$tmp[0]);
            if ($tmp[3]>=($len[-1]*0.8) && $tmp[2]>80 && $tmp[3]>=20){
                if ($tmp[8]>$tmp[9]){
                    $begin=$tmp[9];
                    $end=$tmp[8];
                }else{
                    $begin=$tmp[8];
                    $end=$tmp[9];
                }
                $end=$end+($len[-1]-$tmp[7]);
                $begin=$begin-($tmp[6]-1);
                if (exists $crispr_hit{$tmp[1]}){
                    $crispr_hit{$tmp[1]}=$crispr_hit{$tmp[1]}.":".$begin."-".$end;
                }else{
                    $crispr_hit{$tmp[1]}=$begin."-".$end;
                }
            }
        }
        close(FILE);
    }
}
closedir(DIR);

open (FILE,"$crisprone/local/bact-draft-repeat-hit.txt");
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if (exists $crispr_hit{$tmp[0]}){
        $crispr_hit{$tmp[0]}=$crispr_hit{$tmp[0]}.":".$tmp[1];
    }else{
        $crispr_hit{$tmp[0]}=$tmp[1];
    }
}
close(FILE);

open (PAL,"$crisprone/local/bacteria-plasmid.list");
%bact_plasmid=();
while (<PAL>){
    @tmp=split(/\s+/,$_);
    if ($tmp[1]=~/plasmid/){
	$bact_plasmid{$tmp[0]}=1;
    }
}
close(PAL);

$new=$ARGV[0]."/spacer-mge.result.new";
open (FILE,"$ARGV[0]/spacer-mge.result");
open (NEW,">$new");
%hit_hash=();
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if ($tmp[-1]=~/^bacteria-complete$/){
	@bacteria_id_tmp=split(/\|/,$tmp[1]);
	if ($bacteria_id_tmp[3]=~/(.+)\.\d+$/){
	    $short_id=$1;
	    if (exists $bact_plasmid{$short_id}){
		print "$_";
		$_=~s/bacteria-complete/plasmid(ncbi)/;
		print "$_";
	    }
	}
    }
    $record="";
    for ($t=0;$t<=@tmp-4;$t++){
	$record=$record." ".$tmp[$t];
    }
    if (exists $hit_hash{$record}){
	next;
    }else{
	$hit_hash{$record}=1;
    }

    if (exists $crispr_hit{$tmp[1]}){
	$flag=0;
	@each=split(/\:/,$crispr_hit{$tmp[1]});
	for ($i=0;$i<=@each-1;$i++){
	    @pos=split(/\-/,$each[$i]);
	    $dis_min=1000000000;
	    for ($j=0;$j<=@pos-1;$j++){
		for ($col=8;$col<=9;$col++){
		    $dis_min_tmp=abs($tmp[$col]-$pos[$j]);
		    if ($dis_min_tmp<$dis_min){
			$dis_min=$dis_min_tmp;
			if ($dis_min<=10){
			    $flag=1;
			}
		    }
		}
	    }
	    if ($flag==1){
		next;
	    }
	}
	if ($flag==0){
	    print NEW "$_";
	}
    }else{
	print NEW "$_";
    }
}
close(FILE);
close(NEW);
`rm $new_repeat`;
`rm $ARGV[0]/*.repeat.out`;
`mv $ARGV[0]/spacer-mge.result $ARGV[0]/spacer-mge.result.old`;
`mv $new $ARGV[0]/spacer-mge.result`;
`sort -nr -k 12 $ARGV[0]/spacer-mge.result > $ARGV[0]/spacer-mge.result.sort`;
`mv $ARGV[0]/spacer-mge.result.sort $ARGV[0]/spacer-mge.result`;
