#......................Packages Used..........................
use Getopt::Std;
use Cwd;
use strict;
use Bio::DB::GenBank;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Basename;


open FH1,">",$ARGV[1];
print FH1 "";
close FH1;
open FH1,">>",$ARGV[1];
open FH2,">",$ARGV[4];
print FH2 "";
close FH2;
open FH2,">>",$ARGV[4];

my $blast_db=$ARGV[2];
my $str;
my $blastout;
$str = Bio::SeqIO->new(-file=>"$ARGV[0]" , -format => 'fasta' );
while (my $input = $str->next_seq()){
	my $identity=0;
	my $QueryID=$input->id;
	my $QuerySeq=$input->seq();
	$QuerySeq=~s/\.//g;
	$QuerySeq=~s/\-//g;
	my $QueryLength=length($QuerySeq);
	#print FH1 $QueryID."\n";
	#print FH1 "$QueryID\t";
	#......................
 my $Blastquery_fasta = basename($ARGV[0])."_VirDB_query.fasta";
	open Blastquery, ">", $Blastquery_fasta;
	print Blastquery ">".$QueryID."\n".$QuerySeq;
	close Blastquery;
	my $thread=int($ARGV[3]);
	$blastout=$Blastquery_fasta.".blastout";
	if (-e $blastout){unlink($blastout);}
 #	system "del blastout";
	#print  "$ARGV[5] -db $blast_db -query $Blastquery_fasta -out $blastout -num_threads $thread\n\n";
	system "$ARGV[5] -db $blast_db -query $Blastquery_fasta -out $blastout -num_threads $thread";#-v 1 -b 1 #-y 60 -Z 150 -W 20 -X 50 -A 10000000";#-W 7 -e 1000 -q -1 -G 1 -E 2 -I T";
	my $par_a=0;
	my $monitor=0;
 	my $qflage=1;
	my $report_obj = new Bio::SearchIO(-format => 'blast',-file   => $blastout);
	while(my $result = $report_obj->next_result ){
		while( my $hit = $result->next_hit ){
			my $hit_acc=$hit->accession;
			my $desc=$hit->description();
		 my $hit_length=$hit->length();
			while( my $hsp = $hit->next_hsp ){
				my $align=$hsp->get_aln();
				my $query_seq="";
				my $Hit_Seq="";
				my $hsp_length=$hsp->length('hit');
				my $idn=$hsp->percent_identity;
				if($idn>= 80 && $hsp_length>=($QueryLength*0.80)&& $hsp_length>= 40&&$qflage==1){
					foreach my $seq ($align->each_seq) {
						if($query_seq){$Hit_Seq= $seq->seq;}
						else{$query_seq=" ";
						 $query_seq=$seq->seq;}
					}
					my $percent_length=$hsp_length/$hit_length*100;
					print  FH2 $QueryID."\t".$desc."\n";
					print  FH1 $QueryID."\t".$hit_acc."\t".int($percent_length)."%\t".$hit_length."\t".int($idn)."%\t".$hsp_length."\t".$desc."\t".$query_seq."\t".$Hit_Seq."\n";
				 	$qflage=0;
				}
			}
		}
	}
 #	print FH1 "\n";
}
if (-e $blastout){unlink($blastout);}

close FH1;
close FH2;
