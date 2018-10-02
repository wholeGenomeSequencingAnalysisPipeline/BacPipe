#!/usr/bin/perl
use strict;
use warnings;
use Excel::Writer::XLSX;
print $ARGV[0]."\t".$ARGV[1]."\t".$ARGV[2]."habaaal\n";
$ARGV[2]="0,1,2,3,4,5,6,7";
my $out=$ARGV[0];
my $ID=$ARGV[1];

my $prokka_gnbnk=$out."/".$ID."/genome_annotation/".$ID.".gbk";

my $prokka=$out."/".$ID."/genome_annotation/results.txt";
my $mlst=$out."/".$ID."/mlst_typing/results.txt";
my $plasmids=$out."/".$ID."/plasmids/results.txt";
my $resistance_profile=$out."/".$ID."/resistance_profile/results.txt";
my $quast=$out."/".$ID."/spades_assembly/quast/report.txt";
my $virulence_profile=$out."/".$ID."/virulence_profile/results.txt";
my $ResFam=$out."/".$ID."/Resfams_directory/results.txt";
my $VirDB=$out."/".$ID."/VirDBSearch_directory/results.txt";
my $card=$out."/".$ID."/CARDsearch_directory/results.txt";
#my $card=$out."/".$ID."/parSNP_directory/results.txt
mkdir($out."/Summary");


my @tools_name= ("mlst","plasmids","resistance_genes","virulence_genes","prokka","ResFam","card","VirDB");
my @tools=($mlst,$plasmids,$resistance_profile,$virulence_profile,$prokka,$ResFam,$card,$VirDB);

my $tmp=$ARGV[2];

my @arr=split(",",$tmp);
my @tools_select;
my @tools_names_select;
for(my $i=0;$i<scalar(@arr);$i++){
	push(@tools_select,$tools[$arr[$i]]);
	push(@tools_names_select,$tools_name[$arr[$i]]);
}


my $excl=$out."/Summary/".$ID.".xlsx";
unlink($excl);
my $workbook  = Excel::Writer::XLSX->new("$excl");

for(my $i=0;$i<scalar(@tools_select);$i++){
   

    open FH,$tools_select[$i];
    my @results=<FH>;
    close FH;
  
    my $worksheet = $workbook->add_worksheet($tools_names_select[$i]);
    
    for(my$j=0;$j<scalar(@results);$j++){
        chomp($results[$j]);
        my $tmp=$j+1;
        my $A="A".$tmp;
        $worksheet->write($A,$results[$j]);
    }

    
}

$workbook->close;



