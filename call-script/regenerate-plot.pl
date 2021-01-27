#!/usr/local/perl

#0: infor file
#1: plot folder

%false=();
open (FILE,"$ARGV[0]");
while (<FILE>){
    if ($_=~/^\#(.+)/){
	$infor=$1;
	@tmp=split(/\s+/,$1);
	if ($tmp[0]=~/\:repeat_pos\:/){
	    if ($tmp[1]=~/tandem-repeat/ || $tmp[2]=~/\(div\-\)/ || $tmp[3]=~/mock/){
		$false{$tmp[0]}=1;
	    }
	}elsif ($tmp[0]=~/\:cas_pos\:/){
	    if ($tmp[0]=~/\cas12/) { }
	    elsif ($tmp[0]=~/\cas13/) { }
	    elsif ($tmp[1]==0 && $tmp[2]<3){
		$false{$tmp[0]}=1;
	    }
	}
    }
}
close(FILE);

opendir (DIR,"$ARGV[1]");
while ($file=readdir(DIR)){
    if ($file=~/(.+)\.plot.txt$/){
	$name=$1;
	$new_file=$ARGV[1]."/".$file.".new";
	$suspicious_file=$ARGV[1]."/".$name.".suspicious";
	open (FILE,"$ARGV[1]/$file");
	open (NEW,">$new_file");
	open (SUS,">$suspicious_file");
	while (<FILE>){
	    @tmp=split(/\s+/,$_);
	    if ($tmp[4]=~/^repeat\:/){ #gi|116626972|ref|NC_008532.1|_1:repeat_pos:649125-650217;gi|116626972|ref|NC_008532.1|_1:cas_pos:COG0640:csa3:1495864-1496649
		$name_tag="repeat_pos";
	    }elsif ($tmp[4]=~/antiRepeat/){
	    }elsif ($tmp[4]=~/^unk$/){
		$name_tag="cas_pos:UNK:UNK";
	    }
	    else{
		$name_tag="cas_pos:".$tmp[4];
	    }
	    $full_name=$tmp[0].":".$name_tag.":".$tmp[1]."-".$tmp[2];
	    #print "$full_name\n";
	    if (!exists $false{$full_name}){
		print NEW "$_";
		#print "$_";
	    }else{
		print SUS "$_";
	    }
	}
	close(NEW);
	close(SUS);
	close(FILE);
	$old_file=$file;
	$old_file=~s/.plot.txt/.plot.ori/g;
	`mv $ARGV[1]/"$file" $ARGV[1]/"$old_file"`;
	`mv $ARGV[1]/"$file.new" $ARGV[1]/"$file"`;
	if(not (-s "$ARGV[1]/$file")) { `rm -f $ARGV[1]/$file`; } #remove empty *.plot.txt file, Ye Nov 2016
    }
}
closedir(DIR);
