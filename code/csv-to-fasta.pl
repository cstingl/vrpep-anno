#!/usr/bin/perl
# MLT/CSV to FASTA (christoph@qfwfq.eu, 2018)

my ($csv_file);

my $column_used = "Peptide";
my $sep = ",";

if(!defined $ARGV[0]){
	print "Usage: csv-to-fasta csv-file.txt\n";
	exit;
} elsif (!-e $ARGV[0]){
	$ARGV[0] =~ s/\\/\//g;
	print "csv-file.txt ".$ARGV[0]." not found or does not exist!\n";
	exit;
} else {
	$ARGV[0] =~ s/\\/\//g;
	$csv_file = $ARGV[0];
	print "> csv file: $csv_file\n";
}

my $out_file = $csv_file;
   $out_file =~ s/\.csv$/\.fasta/i;
   $out_file =~ s/\s/_/g;

my $INQ_file = $csv_file;
   $INQ_file =~ s/\.csv$/_INQ\.fasta/i;
   $INQ_file =~ s/\s/_/g;


my $INQ_idx  = $csv_file;
   $INQ_idx  =~ s/\.csv$/_INQ\.txt/i;
   $INQ_idx  =~ s/\s/_/g;


if($out_file eq $csv_file){
	print "!! Input file cannot be a fasta file ()csv expected!!\n";
	exit;
}

my %entries = ();
my $peptides=();
my $c = -1;

open(CSV, $csv_file) || die $!;
open(OUT, ">".$out_file) || die $!;
open(INQ, ">".$INQ_file) || die $!;
open(IDX, ">".$INQ_idx ) || die $!;

print IDX join("\t", ("accnr", "pepseq", "n_var", "n_occ"))."\n";

my $header = <CSV>;
chomp($header);
my @H = split(/$sep/, $header);
foreach my $i (0..$#H){
	$c = $i if $H[$i] eq $column_used;
}

if($c == -1){
	print "!! Invalid format (no column '$column_used' found.) !!\n";
	print " -- ".join("\n -- ", @H)."\n\n";
	exit;
} else {
	print "'$column_used' found in column $c.\n";
}

my $show_examples = 10;
while (my $l = <CSV>){
	my @L = split(/$sep/, $l);
	my ($prec, $pep_mod, $pep_INQ );
	if ($L[$c] =~ />/){
		($prec, $pep_mod) = split(/>/, $L[$c]);
	} else {
		$pep_mod  = $L[$c];
	}

	$INQ_flg = 0;
	my $pep_INQ = $pep_mod;
	if($pep_INQ =~ /\(/){
		if($pep_INQ =~ /N\(/ || $pep_INQ =~ /Q\(/ ){
			$pep_INQ =~ s/N\(\+\.98\)/D/g;
			$pep_INQ =~ s/Q\(\+\.98\)/E/g;
			$INQ_flg = 1;
		} 
		if($pep_INQ =~ /I/){
			$pep_INQ =~ s/I/L/g;
			print "!! isoleucin in peptide $pep_INQ found! Replacing by Leucing (peaks denovo default.)\n";
			<STDIN>;
			$INQ_flg = 1;
		}
		$pep_INQ =~ s/[0-9+\-#.\[\]\(\)]//g;
	}

	my $pep = $pep_mod;
	   $pep =~ s/[0-9+\-#.\[\]\(\)]//g;

	if($show_examples > 0 && $pep ne $pep_mod){
		print " > removing extra symbols: $pep -> ";
		$pep =~ s/[0-9+\-#.\[\]\(\)]//g;
		print "$pep (just first $show_examples shown)";
		if($pep ne $pep_INQ){
			print "~> pep/INQ: $pep_INQ\n";
		} else {
			print "\n";
		}
		$show_examples--;
	}

	$prec = $pep if $prec eq "";

	my $accnr = $prec.".0";
	if (!exists $entries{$accnr}{$pep}){
		$entries{$accnr}{$pep} = 1;
	} else {
		$entries{$accnr}{$pep}++;
	}

	if($INQ_flg  == 1){
		my $accnr = $prec.".X";
		if(!exists $entries{$accnr}{$pep_INQ}){
			$entries{$accnr}{$pep_INQ} = 1;
		} else {
			$entries{$accnr}{$pep_INQ}++;
		}
	}
}

close CSV;

my $max_var = 0;

foreach my $accnr (sort {$a cmp $b} keys(%entries)){
		my $var = 0;
		foreach my $pep (sort {$a cmp $b} keys(%{$entries{$accnr}}))  {
			$var++;
			my $accnr_var = $accnr;
			   $accnr_var =~ s/\.X$/\.$var/;

			if($accnr =~ /\.0$/){
				my $accnr0 =$accnr;
				   $accnr0 =~ s/\.0$//;

				print OUT ">".$accnr0."\n".$pep."\n";
			}

			print INQ ">".$accnr_var."\n";
			print INQ "$pep\n";
			print IDX join("\t", ($accnr_var, $pep, $var, $entries{$accnr}{$pep} ))."\n";
		}
		$var_max = $var if $var > $var_max;
}

close IDQ;
close IDX;
