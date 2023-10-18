#!\bin\perl.exe

my $file = shift(@ARGV);
   $file =~ s/\\/\//g;
print "read from  -> $file\n";

my $out_file = $file;
	$out_file =~ s/\.txt/\.tidy\.txt/i;
print "write to -> $out_file\n";

if($file eq $out_file){
	print "!! textfile expected as input!!\n";
	exit;
}


open(TXT, $file) || die "Cannot open $file.\n";
open(OUT, ">$out_file") || die $!;


#print OUT join("\t", ("query", "source_ID", "alignment", "from", "to", "length", "matches", "mismatches", "gaps", "percent_identity", "query_id","subject_id","perc_identity","alignment_length","mismatches","gap_opens","gaps","q_start","q_end","s_start","s_end","evalue","bit_score"))."\n";
print OUT join("\t", ("query", "alignment", "from", "to", "length", "matches", "mismatches", "gaps", "percent_identity", "subject_id"))."\n";

my @check = ("","","Alignments","");

my ($query);
my $lcnt = 0; my $read_A = 0;
my ($alignment_anno, $qseq_from, $qseq, $qseq_to, $aa_simil, $accnr, $start, $delta_seq, $stop, );
my @alignment = ();

my $entry_cnt = 0;
my $query_cnt = 0;
my $align_cnt = 0;
my $CDR1_cnt = 0;
my $CDR2_cnt = 0;
my $CDR3_cnt = 0;

my $rep_cnt=0;
my $rep_stp=1000;

print "reading ";

while(<TXT>){

	$lcnt++;
	chomp;
	if($_ =~ /^# Query:/){
		$entry_cnt++;
		($query) = $_ =~ /# Query: (.*)/;
		#print ">> QUERY -> $query\n";
	} elsif ($_ =~ /# Alignment/){
		$read_A = 1;
	} elsif ($read_A == 1 && ($_ eq "") ){
		
		$read_A = 0;
	} elsif ($read_A == 1 ){
		$_ =~ s/^V\s+//;
		push(@alignment, $_);
		$CDR1_cnt++ if $_ =~ /CDR1/;
		$CDR2_cnt++ if $_ =~ /CDR2/;
		$CDR3_cnt++ if $_ =~ /CDR3/;
	} elsif ($_ =~ /^#.*hits found$/){
		my $first_hit = <TXT>;
		chomp( $first_hit );
		$first_hit = (split(/\t/,$first_hit))[2];
		$query_cnt++ if scalar(@alignment) > 0;
		foreach my $al (@alignment){
			if($al =~ /N\/A\tN\/A\tN\/A\tN\/A\tN\/A\tN\/A\tN\/A/){
				$al =~ s/N\/A\tN\/A\tN\/A\tN\/A\tN\/A\tN\/A\tN\/A/N\/A\tN\/A\tN\/A\tN\/A\tN\/A/;
				print "!";
			}
			print OUT join("\t", ($query, $al, $first_hit))."\n";
			$align_cnt++;
		}
		$query     = "";
		@alignment = ();
		
		$rep_cnt++;
		if($rep_cnt>$rep_stp){
			print ".";
			$rep_cnt=0;	
		}
		
	}
}

print " $entry_cnt entries, $query_cnt queries alignemt, and $align_cnt alignments (CDR1/2/3 = $CDR1_cnt/$CDR2_cnt/$CDR3_cnt).\n";