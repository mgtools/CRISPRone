#!/usr/local/perl

$result=`java -cp ../bin/metaCRT.jar crt  $ARGV[0] $ARGV[1]`;

@tmp=split(/\n/,$result);

$check_result="";
for ($i=0;$i<=@tmp-1;$i++){
    if ($tmp[$i]=~/Found/){
	$check_result=$check_result."\n".$tmp[$i];
    }
}

if ($check_result=~/Found/){
    $check_result=$check_result."\nDone";
    print "$check_result\n";
}else{
    print "NO CRISPR array was found!\nDone\n";
}

