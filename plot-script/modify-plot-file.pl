#!/usr/local/perl

opendir (DIR,"$ARGV[0]");
while ($file=readdir(DIR)){
    if ($file=~/(.+).plot.txt/){
	$id=$1;
	open (FILE,"$ARGV[0]/$file");
	@line=<FILE>;
	close(FILE);
	$new=$ARGV[0]."/".$file.".new";
	@tmp=split(/\s+/,$line[0]);
	$start=$tmp[1];
	$end=$tmp[2];
	open (NEW,">$new");
	#print NEW "$line[0]";
	print NEW "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\n";
	if (@line>1){
	    for ($i=1;$i<=@line-1;$i++){
		@tmp=split(/\s+/,$line[$i]);
		if ($tmp[2]>$start && $tmp[2]<$end){
		}else{
		    $start=$tmp[1];
		    $end=$tmp[2]; 
		    print NEW "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\n";
		}
		
	    }
	}
	close(NEW);
	`mv "$new" $ARGV[0]/"$file"`;
    }
}
closedir(DIR);
