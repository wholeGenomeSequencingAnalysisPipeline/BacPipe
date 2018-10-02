open FH,"$ARGV[0]";
my @faa=<FH>;
close FH;
my $mark=0;
unlink("hypothetical_protein.faa");

open FHH,">>hypothetical_protein.faa";

for(my $i=0;$i<scalar(@faa);$i++){
    if($faa[$i]=~/hypothetical protein/){
        $mark=1;
        print FHH $faa[$i];
    }else{
        if($mark==1){
            if($faa[$i]=~/>/){$mark=0;}
            else{print FHH $faa[$i];}
        }
    }
}
close FH;
close FHH;

system "/Users/mohamedahmed/Desktop/UA/WGS_pipeline/hmmer-3.1b2-macosx-intel/binaries/hmmscan --cut_ga --tblout results Resfams.hmm ./hypothetical_protein.faa";

#unlink $ARGV[1];
#open FH,">>$ARGV[1]";
