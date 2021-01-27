#!/usr/local/perl

opendir (DIR,"$ARGV[0]");
while ($pdf=readdir(DIR)){
    if ($pdf=~/(.+)\.pdf/){
	$png=$ARGV[0]."/".$1.".png";
	`convert -trim -density 150 -quality 100 "$ARGV[0]/$pdf" "$png"`;
    }
}
closedir(DIR);
