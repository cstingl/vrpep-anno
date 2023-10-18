#!/usr/bin/perl -w
use strict;
$|=1;

my $tag_min_len = 3;
my %columns = (file=>1, scan=>4, modpepseq=>3, localconf=>17, denovoscore=>6, ALC=>7);

foreach (0..$#ARGV){$ARGV[$_] =~ s/\\/\//g};


unless(defined $ARGV[0] && -e $ARGV[0]  && $ARGV[0] =~ /de.novo.peptides.csv$/i ){
    print "Usage: peaks-anno-localconf.pl \"de novo peptides.csv\" (local conf range)\n";
    exit;
}

my $denovo = shift(@ARGV);
my @lconf_range = (80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99);

my $out_file = $denovo;
   $out_file =~ s/\.csv/\.TAGS\.txt/g;
   $out_file =~ s/de.novo.peptides/de-novo-peptides/g;

my $fasta_file = $denovo;
   $fasta_file =~ s/\.csv/-TAGS\.fasta/g;
   $fasta_file =~ s/de.novo.peptides/de-novo-peptides/g;

print "< reading from $denovo & \n";
open(DENOVO, $denovo) || die "!! Cannot read the denovo  file '$denovo' !!\n";


print "\n";
my $header = <DENOVO>;
chomp($header);
my @H = split(/\,/, $header);
print $header."\n";
foreach (0..$#H){
	print "$_ -> ".$H[$_]."\n" ;
}


print "> writing to $fasta_file.\n ";
print "> writing to $out_file ";
open(OUT, ">$out_file") || die "!! Cannot write to $out_file !!\n";
print OUT join("\t", ("file", "scan", "modpepseq", "seqlen", "denovoscore", "ALC", "localconf", "tag", "tag_start", "tag_end", "tag_len"))."\n";

open(FASTA, ">$fasta_file") || die "!! Cannot write to $fasta_file !!\n";

my $ecnt = 0;
my $scnt = 0;
my $start = time();

my %fasta = ();

while(my $line = <DENOVO>) {
	chomp($line);
	my @L = split(/\,/, $line);
	my ($file) = $L[$columns{file}] =~ /(.*)\..*/;
	my ($scan) = $L[$columns{scan}] =~ /.*\:(.*)/;
	my ($modpepseq) = $L[$columns{modpepseq}];

	$modpepseq =~ s/N\(\+\.98\)/n/g;
	$modpepseq =~ s/Q\(\+\.98\)/q/g;
	$modpepseq =~ s/C\(\+45\.99\)/c/g;
	$modpepseq =~ s/M\(\+15\.99\)/m/g;

	if($modpepseq =~ /\(/){
			print "> UNEXPECTED MODIFICATION IN $modpepseq\n"; 
			exit;
	}

	my @PEPSEQ = split(//, $modpepseq);

	my @localconf = split(/\s+/, $L[$columns{localconf}]);

	if ($#PEPSEQ != $#localconf){
		print "PEPSEQ and LOCALCONF: DIFF LENGTH ($#PEPSEQ != $#localconf) !!!\n";
		print "  $modpepseq (-> ".join(" ", @PEPSEQ).") ~ ".$L[$columns{localconf}]." (".join("|", @localconf).")\n";
		exit;
	}
	foreach my $conf_thresh (@lconf_range){
		my $read = 0;
		my @TAGS = ();
		my @TAG1 = ();
		my @TAGn = ();

		foreach my $pos (0..$#PEPSEQ){
	#		print "[$pos]: ".$PEPSEQ[$pos]." = ".$localconf[$pos]."% -> ";
			if($read == 0){
				if($localconf[$pos] >= $conf_thresh ){
					push(@TAGS, $PEPSEQ[$pos]);
					push(@TAG1, $pos+1);
					$read = 1;
				} 
			}  else {
				if($localconf[$pos] >= $conf_thresh ){
					$TAGS[$#TAGS] .= $PEPSEQ[$pos];
				} else {
					push(@TAGn, $pos);
					$read = 0;
				}

			}
			push(@TAGn, $pos+1) if $read == 1 && $pos == $#PEPSEQ ;
		}

		foreach my $i (0..$#TAGS){
			my $tag_len = length($TAGS[$i]);
			my $seq_len = length($modpepseq);
			if($tag_len >= $tag_min_len){
				print OUT join("\t", ($file, $scan, $modpepseq, $seq_len, $L[$columns{denovoscore}], $L[$columns{ALC}], $conf_thresh, $TAGS[$i], $TAG1[$i], $TAGn[$i]), $tag_len )."\n";

				my $curr_tag = uc($TAGS[$i]);
				unless(defined $fasta{$curr_tag}){
					print FASTA ">".$curr_tag."\n".$curr_tag."\n";
					$fasta{$curr_tag} = 1;
				} else {
					$fasta{$curr_tag}++;
				}
			}
		}
	}

	$ecnt++;
	$scnt++;
	if($scnt>=1000){
		print ".";
		$scnt=0;
	}

}

my $run_time = time() - $start;

print " done (in $run_time sec.).\n";

close DENOVO;
close FASTA;
