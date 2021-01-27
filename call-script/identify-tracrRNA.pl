#!/usr/local/perl

#0: plot folder
#1: repeat file
#2: fasta file



%cas_loci=();
%cas9_dir=();
%all_cas=();
opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/(.+)\.plot\.txt$/){
	$sample=$1;
	open (FILE,"$ARGV[0]/$file");
	@line=<FILE>;
	close(FILE);
	for ($i=0;$i<=@line-1;$i++){
	    chomp $line[$i];
	    @cas_element=split(/\s+/,$line[$i]);
	    if (exists $all_cas{$cas_element[0]}){
		$all_cas{$cas_element[0]}=$all_cas{$cas_element[0]}.";".$line[$i];
	    }else{
		$all_cas{$cas_element[0]}=$line[$i];
	    }

	    if ($line[$i]=~/\s+repeat\:(.+)\s+r\d+/){
		$repeat_len=length($1);
		@repeat=split(/\s+/,$line[$i]);
		
		$system_start=$repeat[1];
		$system_end=$repeat[2];
		$flag=0;
		$dis_flag=0;
		$key=$sample.":repeat_pos:".$repeat[1]."-".$repeat[2].":".$repeat_len;
		$check_cas9=0;
		
		if ($i!=0){
		    $pre=$i-1;

		    @pre_tmp=split(/\s+/,$line[$pre]);
		    $max_pre_dis=$repeat[1]-10000;
		    if($max_pre_dis>0 && $pre_tmp[2]<$max_pre_dis){
			$system_start=$repeat[1];
			$dis_flag=1;
			#$flag=1;
			
		    }else{
		    if($pre_tmp[2]>=($repeat[1]-1500)){ 
			if ($i<=2){
			    for ($j=($i-1);$j>=0;$j--){
				if (($line[$j]=~/type\-II\-/) || ($line[$j]=~/Type\-II$/)){
				    $flag=1;
				    $pre=$j;
				    last;
				}
			    }
			}else{
			    $count=0;
			    $j=$i-1;
			    while ($count<=2){
				if ($line[$j]=~/unknown$/){
				    $j--;
				    next;
				}else{
				    $count++;
				    if ($line[$j]=~/type\-II\-/ || $line[$j]=~/Type\-II$/){
					
					$flag=1;
					$pre=$j;
					last;
				    }
				    $j--;
				}
				
			    }
			}
			if ($flag==1){
			    @pre_tmp=split(/\s+/,$line[$pre]);
			    if ($pre_tmp[2]>($repeat[1]-10000)){
				$system_start=$pre_tmp[1];
				while (($pre_tmp[5]=~/type\-II\-/ || $pre_tmp[5]=~/Type\-II$/ || $pre_tmp[5]=~/unknown/) && ($pre<$i) && ($pre>=0)){
				    
				    if ($pre_tmp[4]=~/\:cas9$/){
					$cas9_dir{$key}=$pre_tmp[3];
					$check_cas9=1;
				    }
				    $system_start=$pre_tmp[1];
				    $pre=$pre-1;
				    @pre_tmp=split(/\s+/,$line[$pre]);
				}
			   }			    

			}
		    }
		}
		
		}#end check repeat is in the first line
		if ($i!=(@line-1)){
		    $next=$i+1;
		    @next_tmp=split(/\s+/,$line[$next]);
		    #print "$repeat[2]\t$next_tmp[1]\n";
		    if ($repeat[2]<($next_tmp[1]-10000)){
			#$flag=1;
			$dis_flag=1;
			$system_end=$repeat[2];
			
		    }else{
		    if($repeat[2]>=($next_tmp[1]-1500)){
			if ($i>=(@line-1-3)){
			    for ($j=($i+1);$j<=@line-1;$j++){
				if ($line[$j]=~/type\-II\-/ || $line[$j]=~/Type\-II$/){
				    $flag=1;
                                    $next=$j;
				    last;
				}
			    }
			    
			}else{
			    $count=0;
                            $j=$i+1;
                            while ($count<=2){
				if ($line[$j]=~/unknown$/){
                                    $j++;
				    next;
				}else{
                                    $count++;
				    
                                    if ($line[$j]=~/type\-II\-/ || $line[$j]=~/Type\-II$/){
                                        $flag=1;
                                        $pre=$j;
                                        last;
                                    }
				    $j++;
                                }
                               
                            }

			}
		    }
			
			if ($flag==1){
			    @next_tmp=split(/\s+/,$line[$next]);
			    if ($next_tmp[1]<($repeat[2]+10000)){
			    while (($next_tmp[5]=~/type\-II\-/ || $next_tmp[5]=~/type\-II$/ || $next_tmp[5]=~/unknown/) && ($next>$i) && ($next<=(@line-1))){
				if ($next_tmp[4]=~/\:cas9$/){
				    $cas9_dir{$key}=$next_tmp[3];
				    $check_cas9=1;
				}
				$system_end=$next_tmp[2];
				$next=$next+1;
				@next_tmp=split(/\s+/,$line[$next]);
			    }
			    }
			    
			}
                    }
	
		}#end check repeat is in the last line
		if ($flag==1||$dis_flag==1){
		    
		    $cas_loci{$key}=$system_start." ".$system_end;
		    
		}
		if ($check_cas9==0){
		    $cas9_dir{$key}="no";
		}
		

	    }#end check repeat
	}#end read lines
    }#end check file end with plot.txt
}#end read folder
closedir(DIR);

$size = keys %cas_loci;

if ($size>=1){
    open (CONT,"$ARGV[1]");
    %contig=();
    while (<CONT>){
	chomp $_;
	if ($_=~/\>(.+)/){
	    $head=$1;
	    @tmp=split(/\s+/,$head);
	    $key=$tmp[0].":".$tmp[1];
	    #$contig{$key}="";
	}
	else{
	    chomp $_;
	    $len=length($_);
	    $key=$key.":".$len;
	    $contig{$key}=$_;
	    #$contig{$key}=$contig{$key}.$_;
	}
    }
    close (CONT);

    $output=$ARGV[0]."/typeII.repeat";
    %typeIIrepeat=();
    open (NEW,">$output");
    for $key (keys %cas_loci){
	if (exists $contig{$key}){
	    typeIIrepeat{$key}=$contig{$key};   ###loading type II repeat seq
	    print NEW "\>$key\n$contig{$key}\n";
	    
	}
    }
    close(NEW);
}

$crRNA=$ARGV[0]."/crRNA.repeat";
$blast_out=$crRNA.".out";


$typeII_repeat_size= -s $output;

if ( (-e $output) && ($typeII_repeat_size>0)){

`perl call-script/negative-strain.pl $output > $crRNA`;

`bin/blast+/bin/makeblastdb  -in $ARGV[2]  -dbtype 'nucl' -parse_seqids`;

`bin/blast+/bin/blastn -db $ARGV[2] -query $crRNA -out $blast_out -outfmt 6  -word_size 11`;

####load crRNA for the future sequence comparison
################################################



$crRNA_result=$ARGV[0]."/crRNA.result";
open (NEW,">$crRNA_result");
open (FILE,"$blast_out");
%crRNA_bash=();
while (<FILE>){
    #print "$_";
    @tmp=split(/\s+/,$_);
    if ($cas9_dir{$tmp[0]}=~/no/){
	next;
    }
    
    if ($tmp[0]=~/(.+)\:repeat_pos\:(\d+)\-(\d+)\:(\d+)$/){
	$id=$1;
	$repeat_start=$2;
	$repeat_end=$3;
	
	$repeat_len=$4;
	
	if ($id eq $tmp[1]){
	    #if (($tmp[3]>=10 && ($tmp[3]<$repeat_len*0.6) && $tmp[4]<=2)||(($tmp[3]>=$repeat_len*0.6) && $tmp[4]<=5)){   #parameters for tracrRNA
	    if ($tmp[3]>=$repeat_len*0.6 && $tmp[4]<=5){
		if ($tmp[8]<$tmp[9]){
		    $start=$tmp[8];
		    $end=$tmp[9];
		    		    		    
		    if ($cas9_dir{$tmp[0]}=~/\+/){
			$rna_dir="+";
		    }else{
			$rna_dir="-";
		    }
		}else{
		    $start=$tmp[9];
                    $end=$tmp[8];
		    if ($cas9_dir{$tmp[0]}=~/\+/){
                        $rna_dir="-";
                    }else{
                        $rna_dir="+";
                    }

		}
		@repeat_array=split(/\s+/,$cas_loci{$tmp[0]});
		if (($start<$repeat_start || $start>$repeat_end) && ($end<$repeat_start || $end>$repeat_end)){
		    
		    if ($start>=($repeat_array[0]-1000) && $start<=($repeat_array[1]+1000) && $end>=($repeat_array[0]-1000) && $end<=($repeat_array[1]+1000)){
			$key=$tmp[1]." ".$start;
			$value=$cas9_dir{$tmp[0]}."\t".$rna_dir."\t".$tmp[1]."\t".$tmp[3]."\t".$tmp[4]."\t".$start."\t".$end;
			$crRNA_bash{$key}=$value;
			print NEW "$cas9_dir{$tmp[0]}\t$rna_dir\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$start\t$end\n";
		    }
		}
	    }
	}
    }
}
close(FILE);
close(NEW);


for $key (keys %all_cas){
    @tmp=split(/\;/,$all_cas{$key});
    my %pos=();
    for ($i=0;$i<=@tmp-1;$i++){
	@tmp_pos=split(/\s+/,$tmp[$i]);
	$pos{$tmp_pos[1]}=$tmp[$i];
    }
    for $cr (keys %crRNA_bash){
	@cr_name=split(/\s+/,$cr);
	if ($cr_name[0] eq $key){
	    @crRNA_order=split(/\s+/,$crRNA_bash{$cr});
	    $pos{$cr_name[1]}=$cr_name[0]."\t".$crRNA_order[5]."\t".$crRNA_order[6]."\t".$crRNA_order[0]."\tantiRepeat\t$crRNA_order[3]-$crRNA_order[4]";
	}
    }
    $output_crRNA=$ARGV[0]."/".$key.".with-crRNA.txt";
    open (NEW,">$output_crRNA");
    for $new ( sort {$a<=>$b} keys %pos) {
	print NEW "$pos{$new}\n";
    }
    close(NEW);
    $old=$ARGV[0]."/".$key.".plot.txt";
    `mv \"$output_crRNA\" \"$old\"`;
}

#`rm $output`;
#`rm $crRNA`;

opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/(.+)\.plot\.txt$/){
	$repeat_all=`grep "repeat:" "$ARGV[0]/$file"`;
	chomp $repeat_all;
	@repeat_tmp=split(/\n/,$repeat_all);
	for ($i=0;$i<=@repeat_tmp-1;$i++){
	    
	}
	$new=$ARGV[0]."/".$file.".clean";
	open (FILE,"$ARGV[0]/$file");
	open (NEW,">$new");
	while (<FILE>){
	    @tmp=split(/\s+/,$_);
	    if ($tmp[4]=~/^antiRepeat$/){
		$flag=0;
		for ($i=0;$i<=@repeat_tmp-1;$i++){
		    @repeat=split(/\s+/,$repeat_tmp[$i]);
		    if (($tmp[1]<$repeat[1] || $tmp[1]>$repeat[2]) && ($tmp[2]<$repeat[1] || $tmp[2]>$repeat[2])){
			
		    }else{
			$flag=1;
			last;
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
	`mv \"$new\" \"$ARGV[0]/$file\"`;
    }
}
closedir(DIR);


#`rm $crRNA_result`;
#`rm $blast_out`;
}
