#!/usr/local/perl

opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/.plot.txt$/){

	$ori=$ARGV[0]."/".$file;
	open(FILE,"$ori");
	$new=$ori."-adjustRNA";
	@line=<FILE>;
	close(FILE);
	@elements=();
	$eleN=0;
	open (NEW,">$new");
	for ($i=0;$i<=@line-1;$i++){
	    @tmp=split(/\s+/,$line[$i]);
	    if ($tmp[4]=~/antiRepeat/){
		$j=$i+1;
		$flag=1;
		@previous=@tmp;
		@next=split(/\s+/,$line[$j]);
		$start=$tmp[1];
		$end=$tmp[2];
		$pre_pos=$i-1;
		if ($pre_pos>=0){
		    @previous_gene=split(/\s+/,$line[$pre_pos]);
		    if ($previous[1]<$previous_gene[2] && $previous[2]<$previous_gene[2]){
			next;
		    }
		}
		while (($next[4]=~/antiRepeat/) && ($next[3] eq $previous[3])){
		    $flag++;
		    $j++;
		    @next=split(/\s+/,$line[$j]);
		    $pre=$j-1;
		    @previous=split(/\s+/,$line[$pre]);
		    $end=$previous[2];
		}
		if ($flag>=4){
		    $elements[$eleN]=$tmp[0]."\t".$start."\t".$end."\t".$tmp[3]."\trepeat:possible\trepeat:".$flag;
		    #print NEW "$tmp[0]\t$start\t$end\t$tmp[3]\trepeat:possible\trepeat:$flag\n";
		    $i=$j;
		    $eleN++
		}else{
		    @currentline=split(/\s+/,$line[$i]);
		    $elements[$eleN]=$line[$i];
		    #print NEW "$line[$i]";
		    $eleN++;
		}
	    }else{
		@currentline=split(/\s+/,$line[$i]);
		$elements[$eleN]=$line[$i];
		#print NEW "$line[$i]";
		$eleN++;
	    }
	}
    
#	close(NEW);
	
#	`rm \"$ori\"`;
#	`mv \"$new\" \"$ori\"`;

	%degrade_repeat=();
	%delete_antirepeat=();
	for ($i=0;$i<=@elements-1;$i++){
	    $pre=$i-1;
	    $after=$i+1;
	    @eachInfor=split(/\s+/,$elements[$i]);
	    
	    if ($eachInfor[4]=~/antiRepeat/){
		
		if ($pre>=0){
		    @eachInfor_pre=split(/\s+/,$elements[$pre]);
		    if($eachInfor_pre[4]=~/^repeat\:/ && abs($eachInfor[1]-$eachInfor_pre[2])<=50){
			$degrade_repeat{$elements[$pre]}="(Degradation:".$eachInfor_pre[2]."bp-".$eachInfor[2]."bp)";
			$delete_antirepeat{$elements[$i]}=1;
			
			next;
		      
		    }
		}
		if ($after<=@elements-1){
		    @eachInfor_after=split(/\s+/,$elements[$after]);
		    if($eachInfor_after[4]=~/repeat\:/ && abs($eachInfor_after[1]-$eachInfor[2])<=50){
			$degrade_repeat{$elements[$after]}="(Degradation:".$eachInfor[1]."bp-".$eachInfor_after[1]."bp)";
                        $delete_antirepeat{$elements[$i]}=1;
                        
			next;
		    }
		}

	    }
	}
	for ($i=0;$i<=@elements-1;$i++){
	    if (exists $degrade_repeat{$elements[$i]}){
		@tmp=split(/\s+/,$elements[$i]);
		print NEW "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]$degrade_repeat{$elements[$i]}\t$tmp[5]\n";
		next;
	    }
	    if (!exists $delete_antirepeat{$elements[$i]}){
		print NEW "$elements[$i]";
	    }
	}
	`mv \"$new\" \"$ori\"`; 

	close(NEW);

    }

}
closedir(DIR);


