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
print FH1 "Gene\t%identity\n";
close FH1;
open FH1,">>",$ARGV[1];
open FH2,">",$ARGV[4];
print FH2 "";
close FH2;
open FH2,">>",$ARGV[4];

system "$ARGV[6] -dbtype nucl -in $ARGV[0]";
my $blast_db=$ARGV[0];
my $str;
my $blastout="./emm_Output.txt";

	my $identity=0;
	my $thread=int($ARGV[3]);
	if (-e $blastout){unlink($blastout);}
	system "$ARGV[5] -p blastn -d $blast_db -i $ARGV[2] -o $blastout -a $thread";#-v 1 -b 1 #-y 60 -Z 150 -W 20 -X 50 -A 10000000";#-W 7 -e 1000 -q -1 -G 1 -E 2 -I T";
	my $par_a=0;
	my $monitor=0;
	#my $qflage=1;
	my $report_obj = new Bio::SearchIO(-format => 'blast',-file   => $blastout);
	while(my $result = $report_obj->next_result ){
		my $QueryID=$result->query_accession;
		my $QueryLength=$result->query_length;
	 my $desc=$result->query_description;
		while( my $hit = $result->next_hit ){
			my $hit_acc=$hit->accession;
		 	my $hit_length=$hit->length();
			while( my $hsp = $hit->next_hsp ){
				my $align=$hsp->get_aln();
				my $query_seq="";
				my $Hit_Seq="";
				my $hsp_length=$hsp->length('hit');
				my $idn=$hsp->percent_identity;
				my $percent_length=$hsp_length/$QueryLength*100;

				if($idn==100 && $hsp_length>=($QueryLength*0.9)){#&& $hsp_length>=($QueryLength*0.97)){
					foreach my $seq ($align->each_seq) {
						if($query_seq){$Hit_Seq= $seq->seq;}
						else{$query_seq=" ";
						$query_seq=$seq->seq;}
					}
					print  FH1 $QueryID."\t".$percent_length."\n";
					print  FH2 $QueryID."\t".$hit_acc."\t".int($percent_length)."%\t".$QueryLength."\t".int($idn)."%\t".$hsp_length."\t".$desc."\t".$query_seq."\t".$Hit_Seq."\n";

				}
			}
		}
	}

if (-e $blastout){unlink($blastout);}
unlink($blast_db.".nhr");
unlink($blast_db.".nin");
unlink($blast_db.".nsq");

close FH1;
close FH2;
