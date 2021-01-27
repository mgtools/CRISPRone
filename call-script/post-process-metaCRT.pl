#!/usr/local/perl


#0: metaCRT result file
open (FILE,"$ARGV[0]");  
@line=<FILE>;
close(FILE);
if ($line[0]=~/^Total\s+sequences\s+checked:/){
    exit;
}


my $folder=$ARGV[0];
if ($folder=~/\//){
    $i=rindex($folder,"/");
    $working_folder=substr($folder,0,$i);
}else{
    $working_folder="./";
}



$new_crt=$ARGV[0].".new";
$old_crt=$ARGV[0];
$log_file=$ARGV[0];
$log_file=~s/.crt/.infor/g;

$check_left_crispr=0;

$crispr_n=1;
for ($i=0;$i<=@line-1;$i++){
    if ($line[$i]=~/^SEQ:\s+(.+)/){
	$head=$1;
	$seq_head=$line[$i].$line[$i+1].$line[$i+2];
    }
    if ($line[$i]=~/^CRISPR\s+\d+\s+Range:\s+(-?\d+)\s+\-\s+(\d+)/){
	$repeat_id=$head."\trepeat_pos:".$1."-".$2;
	my @repeat=();
	my @spacer=();
	$flag=1;
	$j=$i+1;
	$array_n=0;
	#$original_array=$line[$i];
	$original_array=$1." - ".$2."\n";
	while ($flag==1){
	    if ($line[$j]=~/^-?\d+\s+(.+)\s+(.+)\s+\[\s+\d+,\s+\d+\s+\]$/){
		$repeat[$array_n]=$1;
		$spacer[$array_n]=$2;
		$array_n++;
		
	    }elsif  ($line[$j]=~/^\d+\s+(\w+)/){
		$repeat[$array_n]=$1;
	    }

	    if ($line[$j]=~/^Repeat-con\s+(\w+)/){
		$repeat_consus=$1;
		@repeat_con=split(//,$repeat_consus);
	    }
	    if ($line[$j]=~/^Repeats:\s+\d+\s+Average\s+Length/){
		$end_array=$line[$j];
		$flag=0;
	    }
	    $original_array=$original_array.$line[$j];
	    $j++;
	}
	
	###############distinguish tandem repeat vs. CRISPR arrays
	if (-e "$working_folder/tmp.spacer"){
	    `rm $working_folder/tmp.spacer`;
	}
	for ($r=0;$r<=@spacer-1;$r++){
	#    print "$spacer[$r]\n";
	    `echo ">$r" >> $working_folder/tmp.spacer`;
	    `echo "$spacer[$r]" >> $working_folder/tmp.spacer`;
	    `bin/cd-hit-v4.6.1-2012-08-27/cd-hit -i $working_folder/tmp.spacer -o $working_folder/tmp.spacer.0.7 -c 0.7 -n 5`;
	    $check_n=`grep ">" $working_folder/tmp.spacer.0.7|wc -l`;
	    chomp $check_n;
	    $original_n=@spacer;
	}
	##########################################################
	#if ($original_n>=3){
	    if ($check_n>=5 || $check_n>=($original_n/2)){
		open (NEW,">>$new_crt");
		print NEW "$seq_head\n";
		print NEW "CRISPR $crispr_n Range: $original_array\n";
		close(NEW);
		$crispr_n++;
	    }else{
		next;
	    }
	#}else{
	#    next;
	#}
	
        ###############calculate degenerate scores for both sides###################
	$left_end=0;
	$right_end=0;
	for ($r=0;$r<=@repeat-1;$r++){
	    @repeat_each=split(//,$repeat[$r]);
	    $mutation=0;
	    for($com=0;$com<=@repeat_each-1;$com++){
		if ($repeat_each[$com]!~/^-$/){
		    if ($repeat_each[$com] eq $repeat_con[$com]){
		    }else{
			if ($repeat_each[$com]!~/N/){
			    $mutation++;
			}
		    }
		}
	    }
	   # print "$r\t$mutation\n";
	    $n_array=@repeat;
	    $left_end=$left_end+$mutation*($n_array-$r);
	    $right_end=$right_end+$mutation*($r+1);
	    #print "$left_end\t$right_end\n";
	}
	open (LOG,">>$log_file");
	print LOG "#repeat_infor\t$repeat_id\tleft_score:$left_end\tright_score:$right_end\t";
	if ($left_end<$right_end){
	    print LOG "right\t";
	}elsif($right_end<$left_end) {
	    print LOG "left\t";
	}else{
	    print LOG "undecided\t";
	}
	print LOG "original_#_spacer:$original_n\t70id_filter:$check_n\n";
	close(LOG);
###############################################################################
    }
    if ($line[$i]=~/^Total\s+sequences\s+checked/){
	chomp $line[$i];
	`echo "$line[$i]" >> $new_crt`;
    }elsif ($line[$i]=~/^Total\s+bases\s+checked:/){
	chomp $line[$i];
	`echo "$line[$i]" >> $new_crt`;
    }elsif($line[$i]=~/^Total\s+sequences\s+with\s+predicted\s+CRISPR/){
	chomp $line[$i];
	`echo "$line[$i]" >> $new_crt`;
    }elsif($line[$i]=~/^Time\s+used:/){
	chomp $line[$i];
	`echo "$line[$i]" >> $new_crt`;
    }
    
}

if (-e "$working_folder/tmp.spacer"){
    `rm $working_folder/tmp.spacer*`;
}

$old_crt=~s/crt/metaCRT.raw/g;
`mv "$ARGV[0]" "$old_crt"`;
`mv "$new_crt" "$ARGV[0]"`;

