#!/usr/local/perl
use POSIX qw/floor/; 
use File::Basename;
my $dirs = dirname($0); #July 2017, YY 
my $crisprone="$dirs/../";

$bin="$crisprone/bin"; #Ye Nov 2016
$spacer_tmp=$ARGV[0].".spacer.tmp";
$crispr_region=$ARGV[0].".region";
$crisprcon=$ARGV[0].".cons-repeat";
$log_file=$ARGV[0];
$log_file=~s/.crt/.infor/g;
$remove=1;

open (FILE,"$ARGV[0]"); #.crt file                                                                                           
%each_spacer=();
%crt=(); 
%cons_repeat=();
while (<FILE>){
    if ($_=~/^SEQ:/){
        @tmp=split(/\s+/,$_);
        $sample=$tmp[1];
    }
    @col=split(/\s+/,$_);
    if (@col==2){
        if ($col[0]=~/^\d+/){
	    $repeat=$col[1];
	    $repeat=~s/\-//g;
	    $crt{$key}=$crt{$key}.$repeat;
	}elsif ($_=~/^Repeat-con\s+(.+)/){
	    $cons_repeat{$key}=$1;
	}
	next;
    }
    if ($_=~/CRISPR\s+\d+\s+Range:\s+(.+)\s+\-\s+(.+)/){
	my @spacer=();
	$start=$1;
        $end=$2;
        if ($start<0){
            $start=0;
        }
        $key=$sample.":repeat_pos:".$start."-".$end."";
        next;
    }
    if ($_=~/^\-?\d+\s+(.+)\s+(.+)\s+\[\s+\d+\,\s+\d+\s+\]/){
        $spacer=$2;
        $repeat=$1;
	$repeat=~s/\-//g;
	$spacer=~s/\-//g;
        if (exists $crt{$key}){
	    $crt{$key}=$crt{$key}.$repeat.$spacer;
	    $each_spacer{$key}=$each_spacer{$key}."\n".$spacer;
	}else{
	    $crt{$key}=$repeat.$spacer;
	    $each_spacer{$key}=$spacer;
	}
	
	next;
    }
}
close(FILE);

open (NEW,">$crispr_region");
for $key (keys %crt){
    print NEW ">$key\n";
    my @seq = unpack("(A70)*", $crt{$key});
    for ($i=0;$i<=@seq-1;$i++){
	print NEW "$seq[$i]\n";
    }
}
close(NEW);


open (NEW,">$crisprcon");
for $key (keys %cons_repeat){
    $repeat_len=length($cons_repeat{$key});
    print  NEW ">$key\:$repeat_len\n";
    print  NEW "$cons_repeat{$key}\n";
}
close(NEW);

%spacer_div=();
%strong=(); #Ye, Nov 2016 -- Bacteroides_sp._D2 two classic CRISPRs were annotated as tandem-repeats
for $key (keys %each_spacer){
    #print "$key\n";
    #print "$each_spacer{$key}\n";
    @spacer_tmp=split(/\n/,$each_spacer{$key});
    if (-e $spacer_tmp){
	if($remove==1) { `rm $spacer_tmp`; }
    }
    
    open (NEW,">$spacer_tmp");
    for ($i=0;$i<=@spacer_tmp-1;$i++){
	print NEW ">$i\n";
	print NEW "$spacer_tmp[$i]\n";
    }
    close(NEW);
    `$bin/cd-hit-v4.6.1-2012-08-27/cd-hit -i $spacer_tmp -o $spacer_tmp.0.7 -c 0.7 -n 5`;
    $check_n=`grep ">" $spacer_tmp.0.7|wc -l`;
    chomp $check_n;
    $original_n=@spacer_tmp;
    #spacer_div
    if ($check_n>=5 || $check_n>=($original_n/2)){
	$spacer_div{$key}=$original_n."->".$check_n."(div+)";
	if($check_n>=5) { $strong{$key} = $check_n; } #Ye Nov 2016
#	$spacer_div{$key}="";
    }else{
	$spacer_div{$key}=$original_n."->".$check_n."(div-)";
#	$spacer_div{$key}="div(-)";
    }
}
if($remove == 1) {
 `rm -f $spacer_tmp $spacer_tmp.0.7 $spacer_tmp.0.7.clstr`; #add -f, Ye, Nov 2016
}

@path=split(/\//,$crispr_region);
$file_name=pop @path;
$new_path=join("/",@path);
$tmp=`pwd`;
chomp $tmp;
chdir ($new_path);
`$bin/trf409.legacylinux64 $file_name 2 7 7 80 10 50 500 -f -d -m -h`;
if($remove == 1) {
 `rm -f *.2.7.7.80.10.50.500.mask`; #add -f Ye, Nov 2016
}
chdir($tmp);
$dat=$crispr_region.".2.7.7.80.10.50.500.dat";
%ori_repeat=();
%max_repeat=();
open (FILE,"$dat");
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if ($_=~/^Sequence\:\s+(.+)/){
	$id=$1;
	if ($id=~/repeat_pos\:(\d+)\-(\d+)/){
	    $ori_repeat{$id}=$2-$1;
	    $max_repeat{$id}=0;
	}
    }elsif (@tmp==15 && $_=~/^\d+/){
	if (($tmp[1]-$tmp[0])>$max_repeat{$id}){
	    $max_repeat{$id}=$tmp[1]-$tmp[0];
	}
    }
}
close(FILE);

%tandem=();
for $key (keys %ori_repeat){
    if ($max_repeat{$key}>=floor($ori_repeat{$key}*0.8)){
	$tandem{$key}="tandem-repeat";
    }
}

$mockoutput=$crisprcon.".mock"; 
`$bin/blast+/bin/blastn -db $crisprone/local/mock-CRISPR.fa -query $crisprcon -out $mockoutput -outfmt 6 -word_size 11 -num_threads 8`;

%mock=();
open (FILE,"$mockoutput");
while (<FILE>){
    @tmp=split(/\s+/,$_);
    @id=split(/\:/,$tmp[0]);
    $repeat_len=pop @id;
    $name=join(":",@id);
    if (!exists $mock{$name}){
	if ($tmp[3]>=floor($repeat_len*0.9) && $tmp[4]<=3){
	    $mock{$name}="mock";
	}

    }
}
close(FILE);
if($remove == 1) {
 `rm -f $mockoutput`;
}

open (NEW,">>$log_file");
for $key (keys %cons_repeat){
    print NEW "\#$key\t";
    if((not exists($strong{$key})) and (exists $tandem{$key})){ #Ye, Nov 2016; add strong{$key}
	print NEW "$tandem{$key}\t";
    }else{
	print NEW "-\t";
    }
    print NEW "$spacer_div{$key}\t";

    if (exists $mock{$key}){
	print NEW "$mock{$key}\n";
    }else{
	print NEW "-\n";
    }
}
close(NEW);
if($remove == 1) {
  `rm -f $crispr_region $crisprcon`;
}
