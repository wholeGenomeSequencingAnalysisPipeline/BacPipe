use strict;
#use Data::Dump qw(dump);
open FHH,">",$ARGV[1];print FHH "";close FHH;
open FHH,">>",$ARGV[1];
my @arr=split(',',$ARGV[0]);
my @matrix;
$matrix[0][0]="Gene";
#@matrix=( [ "one", "two", "three" ],[  4,   5,  6],['habal']);
my $mlst="";

if ($ARGV[2] eq "card"){
for(my $i=0;$i<scalar(@arr);$i++){
	my @habal=split("/",$arr[$i]);
 
	$matrix[0][$i+1]=$habal[scalar(@habal)-3];

	#print $arr[$i]."\n";
	open FH,$arr[$i];
	my @file=<FH>;
	close FH;
	for(my $j=0;$j<scalar(@file);$j++){
	 	if($file[$j]=~/^#/){$j=scalar(@file);next;}
		chomp($file[$j]);
		my @tmp=split("	",$file[$j]);
		my $marker=0;
	 	#print $#matrix."\n";
		for(my $o=0;$o<$#matrix+1;$o++){
		 #print $tmp[0]."	".$matrix[$o][0]."\n";
			if($tmp[1] eq $matrix[$o][0]){
			 	$marker=1;
			 	my $value="Y";#"100%";
				if($matrix[$o][$i+1]){
					$matrix[$o][$i+1]=$matrix[$o][$i+1].",".$value;
				}
				else{
					$matrix[$o][$i+1]=$value;
				}
			}
		}
		if($marker==0){
		 	my $value="Y";#"100%";
			$matrix[$#matrix+1][0]=$tmp[1];
		 	$matrix[$#matrix][$i+1]=$value;
		}
	}
}


}
elsif ($ARGV[2] eq "vfdb"){
for(my $i=0;$i<scalar(@arr);$i++){
	my @habal=split("/",$arr[$i]);
 
	$matrix[0][$i+1]=$habal[scalar(@habal)-3];

	#print $arr[$i]."\n";
	open FH,$arr[$i];
	my @file=<FH>;
	close FH;
	for(my $j=0;$j<scalar(@file);$j++){
	 	if($file[$j]=~/^#/){$j=scalar(@file);next;}
		chomp($file[$j]);
		my @tmp=split("	",$file[$j]);
	 	my @tmp_=split('\)',$tmp[1]);
	 	$tmp_[0]=~s/^\(//g;
		my $marker=0;
	 	#print $#matrix."\n";
		for(my $o=0;$o<$#matrix+1;$o++){
		 #print $tmp[0]."	".$matrix[$o][0]."\n";
			if($tmp_[0] eq $matrix[$o][0]){
			 	$marker=1;
			 	my $value="Y";#"100%";
				if($matrix[$o][$i+1]){
					$matrix[$o][$i+1]=$matrix[$o][$i+1].",".$value;
				}
				else{
					$matrix[$o][$i+1]=$value;
				}
			}
		}
		if($marker==0){
		 	my $value="Y";#"100%";
			$matrix[$#matrix+1][0]=$tmp_[0];
		 	$matrix[$#matrix][$i+1]=$value;
		}
	}
}
 
 
}
else{
for(my $i=0;$i<scalar(@arr);$i++){
 my @habal=split("/",$arr[$i]);
 
 $matrix[0][$i+1]=$habal[scalar(@habal)-3];

 #print $arr[$i]."\n";
	open FH,$arr[$i];
 	my @file=<FH>;
	close FH;
 	if($file[0]=~/domain number estimation/){
		for(my $j=3;$j<scalar(@file);$j++){
		 	if($file[$j]=~/^#/){$j=scalar(@file);next;}
			chomp($file[$j]);
			my @tmp=split(" ",$file[$j]);
			my $marker=0;
		 	#print $#matrix."\n";
			for(my $o=0;$o<$#matrix+1;$o++){
			 #print $tmp[0]."	".$matrix[$o][0]."\n";
				if($tmp[0] eq $matrix[$o][0]){
				 	$marker=1;
				 	my $value="100%";
					if($matrix[$o][$i+1]){
						$matrix[$o][$i+1]=$matrix[$o][$i+1].",".$value;
					}
					else{
						$matrix[$o][$i+1]=$value;
					}
				}
			}
			if($marker==0){
			 	my $value="100%";
				$matrix[$#matrix+1][0]=$tmp[0];
			 	$matrix[$#matrix][$i+1]=$value;
			}
		}
	}
	else{
	 
	if($file[scalar(@file)-1]=~/MLST Profile: ([\w\W]+)$/){
		my $tmp=$1;
		if($file[0]=~/Sequence Type: ([\w\W]+)$/){
		 #print $file[scalar(@file)-1].$file[0];
			my $tmp2=$1;
			chomp($tmp);
			chomp($tmp2);
	 		$matrix[1][$i+1]=$tmp.":".$tmp2;
			if($matrix[1][0]){
				
			}
			else{
				$matrix[1][0]="MLST";
			}
		}
	}
	 my$vir=0;
	 my @tmp=split("\t",$file[0]);
	 if($tmp[0] eq "Database"){$vir=1;}
	 my $x=0;my$y=1;
	 #print $vir."\n";
	 if($vir==1){$x=1;$y=2;}
	for(my $j=1;$j<scalar(@file);$j++){
	 #print $file[$j];
	 chomp($file[$j]);
		my @tmp=split("\t",$file[$j]);
	 #print $x."\t".$y."\t".$tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\n";

	 	if($tmp[0] eq "Gene"){next;}
	 	elsif($tmp[0] eq "Locus"){next;}
	 #elsif($tmp[0] eq "Database"){$vir=1;next;}
	 	elsif($file[$j]=~/Please note /){last;$j=scalar(@file);}
	 	elsif($file[$j] eq ""){last;}
	 	#elsif($tmp[0]=~/Organism:/){last;}
		else{
			#print "!".$tmp[0]."1\n";
		 	my $marker=0;
		 	#print $#matrix."\n";
		 #my $x=0;my$y=1;
		 #print $vir."\n";
		 #if($vir==1){$x=1;$y=2;}
			for(my $o=0;$o<$#matrix+1;$o++){
			 
			 #print $tmp[0]."	".$matrix[$o][0]."\n";
				if($tmp[$x] eq $matrix[$o][0]){
				 	$marker=1;
				 	my $value=int($tmp[$y])."%";
					if($matrix[$o][$i+1]){
						$matrix[$o][$i+1]=$matrix[$o][$i+1].",".$value;
					}
					else{
					 
						$matrix[$o][$i+1]=$value;
					}
				}
			
			}
			if($marker==0){
			 	my $value=int($tmp[$y])."%";
				$matrix[$#matrix+1][0]=$tmp[$x];
			 	$matrix[$#matrix][$i+1]=$value;
			}
		}
	}
	}
 
}
}

for(my $i=0;$i<$#matrix+1;$i++){
	for(my $j=0;$j<scalar(@arr)+1;$j++){
		if($matrix[$i][$j]){}
		else{$matrix[$i][$j]="-";}
	}
}

for(my $i=0;$i<$#matrix+1;$i++){
 	my $str="";
	for(my $j=0;$j<scalar(@arr)+1;$j++){
		$str=$str.$matrix[$i][$j]."\t";
	}
	chop($str);
 	print FHH $str."\n";
}
close FHH;
#dump(@matrix);
=head
my @arr=( [ "one", "two", "three" ],
[  4,   5,  6,  7  ],
[ "alpha", "beta" ]
);
dump(@arr);
=cut
=head
plasmids for gram negative
Plasmid	Identity	Query/HSP	Contig	Position in contig	Note	Accession no.
Col(BS512)	100.00	193/193	NODE_2_length_2165_cov_5.716	242..434		NC_010656
Col(BS512)	100.00	233/233	NODE_2_length_2165_cov_5.716	506..738		NC_010656



Sequence Type: Unknown ST
Gene	% Identity	HSP Length	Allele Length	Gaps	Best match
adk	91.30	23	536	0	adk-178
fumc	90.48	21	382	0	fumc-798
gyrb	100	33	460	0	gyrb-1
icd	95.24	21	518	0	icd-40
mdh	90.91	22	452	0	mdh-299
pura	91.67	24	478	0	pura-145
reca	100	244	510	0	reca-255
Please note that one or more loci do not match perfectly to any previously registered MLSTallele. We recommend verifying the results by traditional methods for MLST!

Organism: Escherichia coli#1
MLST Profile: ecoli




Resistance gene	Identity	Query/HSP	Contig	Position in contig	Phenotype	Accession no.
aadA1	100.00	792/792	NODE_1762_length_863_cov_188.09	12..803	Aminoglycoside resistance	JN815078
blaTEM-2	100.00	861/861	NODE_23_length_53206_cov_72.6697	13388..14248	Beta-lactam resistance	X54606
cfxA	99.69	966/966	NODE_176_length_10297_cov_10.7416	869..1834	Beta-lactam resistance	U38243
cat	98.63	654/655	NODE_13_length_118190_cov_30.0541	28352..29005	Phenicol resistance	M11587
catA1	99.85	660/660	NODE_23_length_53206_cov_72.6697	46469..47128	Phenicol resistance	V00622
QnrB19	100.00	645/645	NODE_895_length_3075_cov_5.55533	1776..2420	Quinolone resistance	HM146784
tet(J)	98.91	1197/1197	NODE_5_length_217152_cov_39.9748	117351..118547	Tetracycline resistance	AF038993
