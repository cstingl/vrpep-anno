#!/usr/bin/perl -w
use strict;
$| = 1;

my ($tag_file, $localconf, $index, $refresh, $show_stats);

# defaults....
$refresh = 0;
$show_stats = 0;

print "+ QUERY Sequence Tags from de-novo-peptides.TAGS.txt\n";
if(!defined $ARGV[0] || !defined $ARGV[1]){
	print "!! No tag file specified !!\n" if !defined $ARGV[0];
	print "!! No local conf threshold specified !!\n" if !defined $ARGV[1];
	&Usage();
	exit;
} elsif (!-e $ARGV[0]){
	print "!! file does not exists !!\n";
	&Usage();
	exit;

} elsif ($ARGV[1] !~ /\d\d/){
	print "!! localconf '". $ARGV[1]."'is not numeric !!\n";
	&Usage();
	exit;	
} elsif($ARGV[1] >= 100 || $ARGV[1] <= 0 ){
	print "!! localconf '". $ARGV[1]."'is not between 0<localconf<100 !!\n";
	&Usage();
	exit;
} else {
	my @extra = ();
	($tag_file, $localconf, @extra) = @ARGV;

	foreach my $a (@extra){
		$refresh    = 1 if $a =~ /refresh/i || $a =~ /redo/i;
		$show_stats = 1 if $a =~ /show/i || $a =~ /stats/i;
	}
	print " > tags file: $tag_file\n";
	print " > local conf threshold: $localconf\n";

	$index = $tag_file;
	$index =~ s/\.TAGS\.txt/\.INQTAGS$localconf\.index/g;
	if($index eq $tag_file){
		print " ! unexpected tags_file (expected = .TAGS.txt). Exiting ...\n";
	} else {
		print " > index file (out): $index\n";
	}
	print " > parameters: refresh = $refresh, show_stats = $show_stats\n";
}

my %TAG = (); # $INQ -> tag


if(-e $index && $refresh == 0){
	open(IDX, "$index")  || die "Cannot read $index.\n";
	print " ~ READING existing $index ";
	<IDX>;
	my $l_cnt = 0;
	while(<IDX>){
		$l_cnt++;
		chomp;
		my @L = split(/\t/, $_);
		if(!defined $L[4]){
			print "> no INQseq defined in column 4 / line $l_cnt:\n";

			foreach my $i ( 0..$#L ){
				print "!! [$i] ".$L[$i]."\n";
			}
			exit;
		}
		
		my $inqseq = $L[4];
		my $tag = $L[5];
		if(exists $TAG{$inqseq}{$tag}){
			$TAG{$inqseq}{$tag}++;	
		} else {
			$TAG{$inqseq}{$tag} = 1;
		}
			#exit if $l_cnt == 1000;

		print "." if $l_cnt % 10000 == 0;
	}
	print " done ($l_cnt lines).\n";

} else {
	open(TAGS, $tag_file) || die "!! Cannot open/access $tag_file !!";
	open(IDX, ">$index")  || die "Cannot write to $index.\n";
	print IDX join("\t", ("file", "scan", "modpepseq", "pepseq", "tagseq", "INQ"))."\n";

	my $header = <TAGS>;
	chomp($header);
	my @H = split(/\t/, $header);



	my $l_cnt = 0;
	my $e_cnt = 0;


	print " ~ START READING ";
	
	while (<TAGS>){
		chomp;
		my @L = split(/\t/, $_);
		my $modpepseq =$L[2];

		if($L[6] == $localconf){
			my $pepseq = uc($modpepseq);
			my $inqseq  = &convINQ($modpepseq);

			my $tagseq  = uc($L[7]);
			my $inqtag  = &convINQ($L[7]);


			print IDX join("\t", (@L[0..1], $modpepseq, $pepseq, $tagseq, 0))."\n";

			if(exists $TAG{$inqseq}{$tagseq}){
				$TAG{$inqseq}{$tagseq}++;	
			} else {
				$TAG{$inqseq}{$tagseq} = 1;
			}
			$e_cnt++;

			if($pepseq ne $inqseq){
				print IDX join("\t", (@L[0..1], $modpepseq, $inqseq, $inqtag, 1))."\n";

				if(exists $TAG{$inqseq}{$inqtag}){
					$TAG{$inqseq}{$inqtag}++;	
				} else {
					$TAG{$inqseq}{$inqtag} = 1;
				}
				$e_cnt++;
			}
		}
		$l_cnt++;
		print "." if $l_cnt % 50000 == 0;
	}
	print " done ($l_cnt lines).\n";
}


if($show_stats == 1){
	print " > Reporing counts of unique peptide-tag pair counts:\n";
	my $scnt = 0;
	my %stat = ();
	foreach my $inqseq (sort keys(%TAG)){
		foreach my $tag(sort keys(%{$TAG{$inqseq}})){
			my $cnt = $TAG{$inqseq}{$tag};
			if(defined $stat{$cnt}){
				$stat{$cnt}++;
			} else {
				$stat{$cnt} = 1;


			}

		}
	}

	foreach my $c (sort {$a <=> $b} keys(%stat)) {
		if( $c <= 10 ){
			print "\t - $c ".$stat{$c}."\n";
		} elsif ($c == 11 ) {
			print "\t - >10: $c=".$stat{$c};
		} else {
			print " $c=".$stat{$c}." ";
		}
	}
	print "\n";
}




# SUBS - SUBS - SUBS - SUBS - SUBS - SUBS - SUBS - SUBS - SUBS - SUBS - SUBS

sub Usage {
	print "Usage: peaks-query-localconf.pl de-novo-peptides.TAGS.txt localconf\n";
	print "+ de-novo-peptides.TAGS.txt: outpout file from peaks-anno-localconf.pl expected.TAGS\n";
	print "+ localconf: local confidence threshold, numeric and < 100 expected.\n";

}

sub convINQ {
	my $pep_INQ = $_[0];
	if($pep_INQ =~ /n/ || $pep_INQ =~ /q/ ){
		#$pep_INQ =~ s/N\(\+\.98\)/D/g;
		$pep_INQ =~ s/n/D/g;
		$pep_INQ =~ s/q/E/g;
	} 
	if($pep_INQ =~ /I/){
		$pep_INQ =~ s/I/L/g;
		print "!";
	}

	return uc($pep_INQ);

}