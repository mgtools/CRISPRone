#!/usr/local/perl

open (FILE,"$ARGV[0]");#.suspicious file

print "<table style=\"width:900px;text-align:center\">";
print"<tr>";
print "<th>Begin-Position</th>";
print "<th>End-Position</th>";
print "<th>Strand</th>";
print "<th>Annotation</th>";
while (<FILE>){
    @tmp=split(/\s+/,$_);
    print"<tr>";
    print "<td>$tmp[1]</td>";
    print "<td>$tmp[2]</td>";
    if ($tmp[4]=~/^repeat\:(.+)/){
	print "<td>.</td>";
	print "<td>CRISPR array ($tmp[4])</td>";
    }else{
	print "<td>$tmp[3]</td>";
	print "<td>$tmp[4]($tmp[5])</td>";
	
    }
    print "</tr>\n";
}
close(FILE);
print "</table>\n";
