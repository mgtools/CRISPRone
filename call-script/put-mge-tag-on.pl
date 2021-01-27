#!/usr/local/perl

%mge_tag=();

=head
open (FILE,"$ARGV[0]"); #spacer-mge.result                                                                                                                                        
while (<FILE>){
    @tmp=split(/\s+/,$_);
    if ($tmp[0]=~/(.+)\:\d+$/){
        $name=$1;
        if ($name=~/(.+)\_(s\d+)$/){
            $seq=$1;
            $spacer_id=$2;
            if ($tmp[12]=~/phage\(HMP\)$/){
                $ann="phage";
            }elsif ($tmp[12]=~/mge\(HMP\)$/){
                $ann="mge";
            }elsif($tmp[12]=~/plasmid\(/){
                $ann="plasmid";
            }elsif($tmp[12]=~/virus\(host\:/){
                $ann="phage";
            }else{
                $ann="mge";
            }
            if (exists $mge_tag{$seq}){
                $mge_tag{$seq}=$mge_tag{$seq}." ".$spacer_id.":".$ann;
            }else{
                $mge_tag{$seq}=$spacer_id.":".$ann;
            }
        }

    }
}
close(FILE);


opendir (DIR,"$ARGV[1]");
=cut
opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/(.+).plot.txt$/){
	$name=$1;
	$new_file=$ARGV[0]."/".$name.".new.plot";
	
	if (exists $mge_tag{$name}){
	    @tmp=split(/\s+/,$mge_tag{$name});
	    %tmp_hash=();
	    for ($i=0;$i<=@tmp-1;$i++){
		@tag=split(/\:/,$tmp[$i]);
		$tmp_hash{$tag[0]}=$tag[1];
	    }
	    open (FILE,"$ARGV[0]/$file");
	    open (NEW,">$new_file");
	    while (<FILE>){
		@line=split(/\s+/,$_);
		if ($line[4]=~/repeat/){
		    @cri_array=split(/\-/,$line[5]);
		    if($cri_array[0]=~/s\d+/){
			if (exists $tmp_hash{$cri_array[0]}){
			    $array_crsipr=$cri_array[0].":".$tmp_hash{$cri_array[0]};
			}else{
			    $array_crsipr=$cri_array[0].":unk";
			}
		    }else{
			$array_crsipr=$cri_array[0];
		    }

		    for ($j=1;$j<=@cri_array-1;$j++){
			if($cri_array[$j]=~/s\d+/){
			    if (exists $tmp_hash{$cri_array[$j]}){
				$array_crsipr=$array_crsipr."-".$cri_array[$j].":".$tmp_hash{$cri_array[$j]};
			    }else{
				$array_crsipr=$array_crsipr."-".$cri_array[$j].":unk";
			    }
			}else{
			    $array_crsipr=$array_crsipr."-".$cri_array[$j];
			}
		    }
		    print NEW "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$array_crsipr\n";
		}else{
		    print NEW $_;
		}
	    }
	    close(FILE);
	    close(NEW);
	}else{
	    
	    open (FILE,"$ARGV[0]/$file");
            open (NEW,">$new_file");
            while (<FILE>){
                @line=split(/\s+/,$_);
		
		if ($line[4]=~/repeat/){
                    @cri_array=split(/\-/,$line[5]);
                    if($cri_array[0]=~/s\d+/){
			$array_crsipr=$cri_array[0].":unk";
		    }else{
                        $array_crsipr=$cri_array[0];
                    }

                    for ($j=1;$j<=@cri_array-1;$j++){
			if($cri_array[$j]=~/s\d+/){
			    $array_crsipr=$array_crsipr."-".$cri_array[$j].":unk";
			}else{
                            $array_crsipr=$array_crsipr."-".$cri_array[$j];
			}
                    }
                    print NEW "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$array_crsipr\n";
		}else{
                    print NEW $_;
                }
            }
            close(FILE);
            close(NEW);
	}
#	`rm \"$ARGV[1]/$file\"`;
	`mv \"$new_file\" \"$ARGV[0]/$file\"`;
    }
}
closedir(DIR);

