#!/usr/local/perl
use File::Basename;
%cas9=();
%cas1=();
%cas2=();
%csn2=();
%cas4=();
my $dirs = dirname($0); #July 2017, YY 
my $casfile="$dirs/../local/cas-db/cas-CRISPR.GA";
open (FILE, $casfile);
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if ($tmp[2]=~/cas1\:/ || ($tmp[2]=~/cas1\_/)){
	$cas1{$tmp[0]}=0;
#	print "$tmp[0]\n";
    }
    if ($tmp[2]=~/cas2\:/){
	$cas2{$tmp[0]}=0;
    }
    if ($tmp[2]=~/csn2\:/ || $tmp[2]=~/csn2b\:/){
	$csn2{$tmp[0]}=0;
    }
    if ($tmp[2]=~/cas9\:/){
	$cas9{$tmp[0]}=0;
#	print "$tmp[0]\n";
    }
    if ($tmp[2]=~/cas4\:/ || $tmp[2]=~/cas4\_/){
	$cas4{$tmp[0]}=0;
    }

}
close(FILE);

@typeII_cas=["TIGR01865","TIGR03031","TIGR00287","TIGR03639", "PF01867.11","TIGR01866","PF09711.5","TIGR00372","TIGR01573",];
opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/(.+).plot.txt$/){
	$new_output=$ARGV[0]."/".$1.".plot.txt.update";
	open (FILE,"$ARGV[0]/$file");
	@line=<FILE>;
	close(FILE);
	$file_line=@line-1;
	%typeII_pos=();
	
	for ($i=0;$i<=@line-1;$i++){
	    @tmp=split(/\s+/,$line[$i]);
	    @flagA=();
	    #@flagB=();
	    @hmm_id=split(/\:/,$tmp[4]);
	    if (exists $cas9{$hmm_id[0]}){
		#print "$hmm_id[0]\n";
 #if (($tmp[4]=~/TIGR01865/) || ($tmp[4]=~/TIGR03031/)){
		@change_pos=();
		$i_pos=0;
		$change_pos[$i_pos]=$i;
		$i_pos++;
		#if ($tmp[4]=~/TIGR01865/){
		$flagA[0]=1; #cas9
		$flagA[1]=0;
		$flagA[2]=0;
		$flagA[3]=0;
		$flagA[4]=0;
#}else{
		 #   $flagB[0]=1;
		 #   $flagB[1]=0;
		 #   $flagB[2]=0;
		 #   $flagB[3]=0;
		#}
		if ($tmp[3]=~/\+/){
		    $next=$i+1;
		    @next_tmp=split(/\s+/,$line[$next]);
		    
		    while (($next_tmp[4]!~/repeat\:/) && ($next<=$file_line)){
			#print "$next\n";
			if ($next_tmp[4]=~/antiRepeat/){
			    $next=$next+1;
			    @next_tmp=split(/\s+/,$line[$next]);
			    next;
			}
			if ($next_tmp[5]!~/Type-II$/ && $next_tmp[5]!~/type-II-/ && $next_tmp[5]!~/universal/ && $next_tmp[5]!~/other/ && $next_tmp[5]!~/^unknown/ && $next_tmp[5]!~/^putative/ ){
			    last;
			}
			$pre=$next-1;
			@pre_tmp=split(/\s+/,$line[$pre]);
			if ($pre_tmp[2]>=($next_tmp[1]+10000)){
			    last;
			}
			
			if ($next_tmp[3]=~/\-/){
			    $next=$next+1;
			    @next_tmp=split(/\s+/,$line[$next]);
			    next;
			}
			@hmm_id=split(/\:/,$next_tmp[4]);
			if (exists $cas2{$hmm_id[0]}){  #cas2
			#if ($next_tmp[4]=~/\:cas2$/){#cas2
			   # if ($flagA[0]==1){
				$flagA[2]=1;
			    #}
			    #if ($flagB[0]==1){
			#	$flagB[2]=1;
			 #   }
			    $change_pos[$i_pos]=$next;
			    $i_pos++;
			}
			if (exists $cas1{$hmm_id[0]}){
			#if (($next_tmp[4]=~/TIGR00287\:/) || ($next_tmp[4]=~/TIGR03639\:/) || ($next_tmp[4]=~/PF01867.11/)){  #cas1
			    #if ($flagA[0]==1){
			    $flagA[1]=1;
                            #}
                            #if ($flagB[0]==1){
                             #   $flagB[1]=1;
                            #}
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
			}
			if (exists $csn2{$hmm_id[0]}){
			#if (($next_tmp[4]=~/TIGR01866\:/) || ($next_tmp[4]=~/PF09711.5/)|| ($next_tmp[4]=~/cd12217/)){#csn2
			 #   if ($flagA[0]==1){
                                $flagA[3]=1;
			  #  }
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
			}
			if (exists $cas4{$hmm_id[0]}){
			    $flagA[4]=1;
			#if ($next_tmp[4]=~/TIGR00372/){#cas4
			 #   if ($flagB[0]==1){
                           #     $flagB[3]=1;
                          #  }
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
			}
			
			$next=$next+1;
			@next_tmp=split(/\s+/,$line[$next]);
		    }
		    $type="Type-II";
#		    print "$file\n";
#		    for ($f=0;$f<=@flagA-1;$f++){
#			print "$flagA[$f]\t";
#		    }
#		    print "\n";
		    if ($flagA[0]==1){
			if (($flagA[1]==1 && $flagA[2]==1) && $flagA[3]==1 && $flagA[4]==0){
			    $type="Subtype-II-A";
			}
			if (($flagA[1]==1 && $flagA[2]==1) && $flagA[3]==0 && $flagA[4]==0){
			    $type="Subtype-II-C";
			}
			if (($flagA[1]==1 && $flagA[2]==1) && $flagA[4]==1 && $flagA[3]==0){
                            $type="Subtype-II-B";
                        }

		    }
		    #if ($flagB[0]==1){
                     #   if (($flagB[1]==1 || $flagB[2]==1) && $flagB[3]==1){
                      #      $type="subtype-II-B";
                      #  }
                        
                    #}
		    for ($pos=0;$pos<=@change_pos-1;$pos++){
			$typeII_pos{$change_pos[$pos]}=$type;
		    }
		} # +

                if ($tmp[3]=~/\-/){
                    $next=$i-1;
                    @next_tmp=split(/\s+/,$line[$next]);
                    while (($next_tmp[4]!~/repeat/) && ($next>=0)){
			if ($next_tmp[4]=~/antiRepeat/){
                            $next=$next-1;
			    @next_tmp=split(/\s+/,$line[$next]);
                            next;
                        }

                        $pre=$next+1;
                        @pre_tmp=split(/\s+/,$line[$pre]);
                        if ($pre_tmp[1]>=($next_tmp[2]+10000)){
                            last;
                        }
			if ($next_tmp[5]!~/Type-II$/ && $next_tmp[5]!~/type-II-/ && $next_tmp[5]!~/universal/ && $next_tmp[5]!~/other/ && $next_tmp[5]!~/^unknown/ && $next_tmp[5]!~/^putative/){
                            last;
                        }

                        if ($next_tmp[3]=~/\+/){
			    $next=$next-1;
			    @next_tmp=split(/\s+/,$line[$next]);
                            next;
                        }
			@hmm_id=split(/\:/,$next_tmp[4]);
			if (exists $cas2{$hmm_id[0]}){
                        #if ($next_tmp[4]=~/\:cas2$/){
                            #if ($flagA[0]==1){
                                $flagA[2]=1;
                            #}
                            #if ($flagB[0]==1){
                            #    $flagB[2]=1;
                            #}
			    $change_pos[$i_pos]=$next;
                            $i_pos++;

                        }
			if (exists $cas1{$hmm_id[0]}){
                        #if (($next_tmp[4]=~/TIGR00287\:/) || ($next_tmp[4]=~/TIGR03639\:/) || ($next_tmp[4]=~/PF01867.11/)){
                         #   if ($flagA[0]==1){
                                $flagA[1]=1;
                         #   }
                         #   if ($flagB[0]==1){
                         #       $flagB[1]=1;
                         #   }
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
                        }
			if (exists $csn2{$hmm_id[0]}){
                        #if (($next_tmp[4]=~/TIGR01866\:/) || ($next_tmp[4]=~/PF09711.5/) || ($next_tmp[4]=~/cd12217/)){
                         #   if ($flagA[0]==1){
                                $flagA[3]=1;
                          #  }
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
                        }
			if (exists $cas4{$hmm_id[0]}){
                        #if ($next_tmp[4]=~/TIGR00372/){
                         #   if ($flagB[0]==1){
                                $flagA[4]=1;
                          #  }
			    $change_pos[$i_pos]=$next;
                            $i_pos++;
                        }

                        $next=$next-1;
                        @next_tmp=split(/\s+/,$line[$next]);
                    }
                    $type="Type-II";
		    if ($flagA[0]==1){
                        if (($flagA[1]==1 && $flagA[2]==1) && $flagA[3]==1 && $flagA[4]==0){
                            $type="Subtype-II-A";
                        }
                        if (($flagA[1]==1 && $flagA[2]==1) && $flagA[3]==0 && $flagA[4]==0){
                            $type="Subtype-II-C";
                        }
                        if (($flagA[1]==1 && $flagA[2]==1) && $flagA[4]==1 && $flagA[3]==0){
                            $type="Subtype-II-B";
                        }

                    }
                		
		    for ($pos=0;$pos<=@change_pos-1;$pos++){
                        $typeII_pos{$change_pos[$pos]}=$type;
		    }

		} # -

	    } #check cas9 end
	}#all lines

        open (FILE,"$ARGV[0]/$file");
	@line=<FILE>;
        close(FILE);
	open (NEW,">$new_output");
	for ($i=0;$i<=@line-1;$i++){
	    @tmp=split(/\s+/,$line[$i]);
	    if ($tmp[4]!~/repeat\:/){
		#for ($j=0;$j<=@typeII_cas-1;$j++){
		 #   $check=index($tmp[4],$typeII_cas[$j]);
		  #  if ($check!=-1){
		#	print NEW "$line[$i]";
		#	next;
		 #   }
		#}
		if (exists $typeII_pos{$i}){
		    print NEW "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$typeII_pos{$i}\n";
		}else{
		    print NEW "$line[$i]";
		}
	    }else{
		print NEW "$line[$i]";
	    }
	}
	close(NEW);
	`mv \"$new_output\" \"$ARGV[0]/$file\"`;
    } # file check end
}  #read dir end
closedir(DIR);


