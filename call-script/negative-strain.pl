#!/usr/local/perl

open (FILE,"$ARGV[0]");
while (<FILE>){
    if($_=~/>/){
	print "$_";
    }else{
	chomp $_;
	$_=uc($_);
	$change=&comp($_);
	$change=&rev($change);
	print "$change\n";
    }
}
close(FILE);

sub comp{
    my ($str_nu)=@_;

    @nuclo=split(//,$str_nu);

    for ($n=0;$n<=@nuclo-1;$n++){
	if ($nuclo[$n]=~/T/){
            $nuclo[$n]=~ s/T/A/;
	}
        elsif ($nuclo[$n]=~/A/){
            $nuclo[$n]=~ s/A/T/;
        }
        elsif ($nuclo[$n]=~/C/){
            $nuclo[$n]=~ s/C/G/;
        }
        elsif ($nuclo[$n]=~/G/){
            $nuclo[$n]=~ s/G/C/;
        }
    }
    $re_str=join("",@nuclo);
    return ($re_str);
}

sub rev{
    my ($strnu)=@_;
    my @nu=split(//,$strnu);
    @renu=reverse(@nu);
    $new=join("",@renu);
    return($new);
}
