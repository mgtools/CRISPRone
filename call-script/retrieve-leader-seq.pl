#!/usr/local/perl

#0 folder
#1 sample 
#2 leader length
#3 fasta file

$folder=$ARGV[0]."/".$ARGV[1]."/";
$leader_length=$ARGV[2];
$gff=$folder."/".$ARGV[1].".gff";
$seq_infor=$folder."/".$ARGV[1].".infor";
$db_seq=$ARGV[3];
$db_seq_file=$db_seq.".nhr";

open (FILE,"$seq_infor");
%infor=();
%degenerate=();
while (<FILE>){
    if ($_=~/^\#/){
	@tmp=split(/\s+/,$_);
	$key=$tmp[1]."-".$tmp[2];
	$degenerate{$key}=$tmp[5];
    }else{
	@tmp=split(/\s+/,$_);
	$infor{$tmp[0]}=$tmp[1];
    }
}
close(FILE);
opendir(DIR,"$folder/plot/");
while ($file=readdir(DIR)){
    if ($file=~/(.+).plot.txt$/){
	$id=$1;
	$new_leader_file=$folder."/plot/".$id.".leader";
	$log=$folder."/plot/".$id.".leader.infor";
	open (FILE,"$folder/plot/$file");
	@line=<FILE>;
	close(FILE);
	for ($i=0;$i<=@line-1;$i++){
	    @tmp=split(/\s+/,$line[$i]);
	    if ($tmp[4]=~/repeat:/){
		$repeat_key=$tmp[0]."-repeat_pos:".$tmp[1]."-".$tmp[2];
		if (!exists $degenerate{$repeat_key}){
		    next;
		}
		open (LOG,">>$log");
		print LOG "$tmp[0]\t$tmp[1]\t$tmp[2]\t";
		$repeat_begin=$tmp[1];
		$repeat_end=$tmp[2];
		$flag_pre=0;
		$flag_after=0;
		$pre=$i-1;
		@pre_tmp=split(/\s+/,$line[$pre]);
		$after=$i+1;
		@after_tmp=split(/\s+/,$line[$after]);
		
		if ($pre>=0){
		    while ( (($line[$pre]=~/\s+tracrRNA\s+/)||($pre_tmp[3]=~/\-/)) && ($pre>0)){
			$pre=$pre-1;
			@pre_tmp=split(/\s+/,$line[$pre]);
		    }
		    if ($pre==0){
			if ($line[$pre]=~/\s+tracrRNA\s+/ || $pre_tmp[3]=~/\-/){
			    $flag_pre=2;
			}else{
			    $flag_pre=1;
			}
		    }else{
			$flag_pre=1;
		    }
		}else{
		    $flag_pre=2;
		}
		if ($flag_pre==1){
		    $pre_end=$pre_tmp[2]+1;
		}elsif($flag_pre==2){
		    if (($repeat_begin-$leader_length)<0){
                        $pre_end=1;
                    }else{
                        $pre_end=$repeat_begin-$leader_length;
                    }
		}
		if (($repeat_begin-$pre_end)>$leader_length){
		    $pre_end=$repeat_begin-$leader_length;
		}
		$last_line=@line-1;
		if ($after<=$last_line){
		    while ( (($line[$after]=~/\s+tracrRNA\s+/)||($after_tmp[3]=~/\+/)) && ($after<$last_line) ){
                        $after=$after+1;
			@after_tmp=split(/\s+/,$line[$after]);
                    }
                    if ($after==$last_line){
                        if ($line[$after]=~/\s+tracrRNA\s+/ || $after_tmp[3]=~/\+/){
                            $flag_after=2;
                        }else{
                            $flag_after=1;
                        }
                    }else{
                        $flag_after=1;
                    }

		}else{
		    $flag_after=2;
		}
		if($flag_after==1){
		    $after_begin=$after_tmp[1]-1;
		}elsif($flag_after==2){
		    if (($repeat_end+$leader_length)>$infor{$id}){
                        $after_begin=$infor{$id};
                    }else{
                        $after_begin=$repeat_end+$leader_length;
                    }
		}
		if (($after_begin-$repeat_end)>$leader_length){
		    $after_begin=$repeat_end+$leader_length;
		}
		$gff_cmd=`grep -w "$id" $gff`;
		chomp $gff_cmd;
		@gff_result=split(/\n/,$gff_cmd);
		for ($j=0;$j<=@gff_result-1;$j++){
		    @gene=split(/\s+/,$gff_result[$j]);
		    if ($gene[6]=~/\+/){
			if ((($gene[4]+1)>=$pre_end)&&($gene[4]<$repeat_begin)){
			    $pre_end=$gene[4]+1;
			}
		    }
		    if ($gene[6]=~/\-/){
			if ((($gene[3]-1)<=$after_begin)&&($gene[3]>$repeat_end)){
			    $after_begin=$gene[3]-1;
			}
		    }
		}
		$repeat_begin=$repeat_begin-1;
		$repeat_end=$repeat_end+1;
		if ($degenerate{$repeat_key}=~/right/){
		    print LOG "$pre_end-$repeat_begin:leader\tc$repeat_end-$after_begin:nonleader\n";
		    #$left_def="leader";
		    #$right_def="nonleader";
		}elsif ($degenerate{$repeat_key}=~/left/){
		    print LOG "$pre_end-$repeat_begin:nonleader\tc$repeat_end-$after_begin:leader\n";
		    #$left_def="nonleader";
                    #$right_def="leader";

		}else{
		    #$left_def="uncertain";
                    #$right_def="uncertain";
		    print LOG "$pre_end-$repeat_begin:uncertain\tc$repeat_end-$after_begin:uncertain\n";
		}
		if (($repeat_begin-$pre_end)>=35){
		    if (!-e $db_seq_file){
			`bin/blast+/bin/makeblastdb  -in $db_seq  -dbtype 'nucl' -parse_seqids`;
		    }
		    `bin/blast+/bin/blastdbcmd -db $db_seq -dbtype 'nucl' -entry "$id" -range "$pre_end-$repeat_begin" -strand "plus" >> "$new_leader_file"`;
		    #$leader_seq=`bin/blast+/bin/blastdbcmd -db $db_seq -dbtype 'nucl' -entry "$id" -range "$pre_end-$repeat_begin" -strand "plus" `;
		    #@result=split(/\n/,$leader_seq);
		    #print ">$id:$pre_end-$repeat_begin\t$left_def\n";
		    #for ($s=1;$s<=@result-1;$s++){
	#		print "$result[$s]\n";
	#	    }
		}
		if (($after_begin-$repeat_end)>=35){
		    if (!-e $db_seq_file){
			`bin/blast+/bin/makeblastdb  -in $db_seq  -dbtype 'nucl' -parse_seqids`;
		    }
		    `bin/blast+/bin/blastdbcmd -db $db_seq -dbtype 'nucl' -entry "$id" -range "$repeat_end-$after_begin" -strand "minus" >> "$new_leader_file"`;
		    #$leader_seq=`bin/blast+/bin/blastdbcmd -db $db_seq -dbtype 'nucl' -entry "$id" -range "$repeat_end-$after_begin" -strand "minus"`;
		    #print "$id:c$after_begin-$repeat_end\t$right_def\n";
		    #@result=split(/\n/,$leader_seq);
		    #for($s=1;$s<=@result-1;$s++){
	#		print "$result[$s]\n";
         #           }

		}
	    }
	    
	}
	#bin/blast+/bin/blastdbcmd -db $db{$pam[-1]} -dbtype 'nucl' -entry "$pam[1]" -range "$pam[2]" -strand "$pam[3]"
    }
}
closedir(DIR);
