#!/usr/bin/perl -w
use strict;
$| = 1;

my $debug = 0;
my $useINQ = 0;

print "+ annotate Ig regions in BlAST (Christoph Stingl, 18. May 2022)\n";
print "  - regions from IgBLAST against IMGT-V\n";
print "  - same accnrs in BLAST and IgBLAST expected.\n";

foreach ( 0 .. $#ARGV ) { $ARGV[$_] =~ s/\\/\//g; }

# fixed lookup tables: IMGT fasta file and regions (results of IgBLAST of IMGT-V-genes)

my $igblast7_tidy  = "external-data/IMGT-V-seq_v210614.igblast7.tidy.txt";
my $IMGT_fasta_file = "external-data/IMGT-VDJC_v210614.fasta";
my ($blast_result, $INQ_file); 

print "+ Starting ...\n";

#                                   _
#  __ _ _ _ __ _ _  _ _ __  ___ _ _| |_ ___
# / _` | '_/ _` | || | '  \/ -_) ' \  _(_-<
# \__,_|_| \__, |\_,_|_|_|_\___|_||_\__/__/
#          |___/

# 1. input filename (blast file)

if(defined $ARGV[0]){
    $blast_result = shift(@ARGV);
    $blast_result =~ s/\\/\//g;
    unless( $blast_result =~ /_INQ\.blast~IMGT-VDJC_v210614.txt/ || $blast_result =~ /\.blast~IMGT-VDJC_v210614.txt/ || $blast_result =~ /_INQ\.blast~IMGT-VDJC_v210614_I2L.txt/ || $blast_result =~ /\.blast~IMGT-VDJC_v210614_I2L.txt/){
        &usage();
        exit;
    }

    # corresponded INQ file (list of accnr (pepseq.X and INQseq)
    if($blast_result =~ /_INQ\./){
        $INQ_file = $blast_result;
        $INQ_file =~ s/_INQ\.blast~IMGT-VDJC_v210614.*\.txt$/_INQ\.txt/i;
        if($INQ_file eq $blast_result){
            &usage();
            exit;
        } else {
            print "  - INQ index: $INQ_file ";
            if(-e $INQ_file){
                print "..ok (exists).\n"
            } else {
                print " ... cannot be accessed (exiting).\n";
                exit;
            }
        }
        $useINQ = 1;        
    } else {
        $INQ_file = $blast_result;
        $INQ_file =~ s/\.blast~IMGT-VDJC_v210614\.*\.txt$/\.txt/i;
        if(-e $INQ_file){
            $useINQ = 1;            
            print " > Use $INQ_file (from filename)\n";
        } else {
            $useINQ = 0;        
        }
    }
        

} else {
    &usage();
    exit;
}


# 2. query arguments

my $query_mode = 0;
my @query_qseqid = ();
my @query_imgtid = ();

# tags loaded by using parameter localconf=
my %TAG = ();
my $use_tags = 0;
my ($localconf);

if(defined $ARGV[0]){
    foreach my $arg (@ARGV){
        if($arg =~ /query=/){
            ($arg) = $arg =~  /query=(.*)/i;
            my @quey_terms = split(/\,/, $arg);
            #example query=AACCTPKSPSSLSASVGNR.0,|IGKV1D-33*01|,|IGKV1-13*01|
            foreach (@quey_terms){
                if($_=~ /\|/ || /\*/){
                    $_ =~ s/\*/x/;
                    $_ =~ s/\|/ /;
                    push(@query_imgtid, $_);
                } elsif($_=~ /\./) {
                    push(@query_qseqid, $_);
                } else {
                    print " ! Unknown query: * for IMGD-IDs and . for INQ peptide accnr expected !\n";
                    exit;
                }
            }

            if(scalar(@quey_terms) > 0){
                print "  ~ query IMGTID: ".join(", ", @query_imgtid)."\n" if scalar(@query_imgtid) > 0;
                print "  ~ query QSeqID: ".join(", ", @query_qseqid)."\n" if scalar(@query_qseqid) > 0;
                $query_mode = 1;
            }
        } elsif ($arg =~ /localconf=/){
            ($localconf) = $arg =~  /localconf=(.*)/i;
            my $tag_file  = $blast_result;
               $tag_file  =~ s/_INQ\.blast~IMGT-VDJC_v210614.*\.txt$/\.TAGS\.txt/i;

            if(!-e $tag_file){
                $tag_file =~ s/de_novo_peptides/de-novo-peptides/gi;
            }


            my $tag_index = $tag_file;
               $tag_index =~ s/\.TAGS\.txt/\.INQTAGS$localconf\.index/g;

            if(!-e $tag_index && !-e $tag_file){
                print "!! TAGs file $tag_file not found !!\n";
                exit;
            } elsif (!-e $tag_index && -e $tag_file){
                print " > preparing TAG index (using peaks-query-localconf.pl) ...";
                my @log = qx( t:\\testarea\\PEAKS\\peaks-query-localconf.pl $tag_file $localconf );
                if(-e $tag_index){
                    print "done.\n";
                } else {
                    print "failed (see log.txt).\n";
                    open(LOG, "> log.txt") || die $!;
                    print LOG join("\n", @log)."\n";
                    close LOG;
                    exit;
                }
            }

            if(-e $tag_index){
                open(IDX, "$tag_index")  || die "Cannot read $tag_index.\n";
                <IDX>;

                while(<IDX>){
                    chomp;
                    my @L = split(/\t/, $_);
                    my $pepseq = $L[3];
                    my $tagseq = $L[4];
                    if(exists $TAG{$pepseq}{$tagseq}){
                        $TAG{$pepseq}{$tagseq}++; 
                    } else {
                        $TAG{$pepseq}{$tagseq} = 1;
                    }
                }
                close IDX;
            } else {
                print "!! Not tag index present and/or generation failed !!\n";
                exit;
            }

            if(scalar(keys(%TAG)) > 0 ){
                $use_tags = 1;
            } else {
                print "!! TAGS file accessible but not tags imported !!\n";
                exit;
            }

        } elsif ($arg =~ /--debug/){
            $debug = 1;
            print "> STARTING IN DEBUGGING MODE (continue with any enter)";
            <STDIN>;
        } elsif ($arg =~ /IMGTfasta=/i){
            ($IMGT_fasta_file) = $arg =~  /IMGTfasta=(.*)/i;
             $IMGT_fasta_file =~ s/\\/\//g;
            print" > IMGT faste: $IMGT_fasta_file\n";

        } elsif ($arg =~ /vgene=/i){
            ($igblast7_tidy) = $arg =~  /vgene=(.*)/i;
            $igblast7_tidy =~ s/\\/\//g;
            print" > V-gene alignment from $igblast7_tidy\n";
            
        }
    }
}


# 3. Derive output file names from paratmeters and inpt filenames

my $out_igblast7_anno = $blast_result;
   $out_igblast7_anno =~  s/\.txt$/\.anno#LC\.txt/i;
    if($use_tags == 1){
        $out_igblast7_anno =~ s/#LC/$localconf/;
    } else {
        $out_igblast7_anno =~ s/#LC//;   
    }
if($out_igblast7_anno eq $blast_result){
    print " !! ERROR: unexpected input file. .igblast7.tidy.txt expected !!\n";
    exit;
} else {
#    print "> peptide level report (best of blast7.tidy annotated) written to $out_igblast7_anno.\n";
}

my $out_igblast7_vreg = $blast_result;
   $out_igblast7_vreg =~  s/\.txt$/\.vreg#LC\.txt/i;
   if($use_tags == 1){   
        $out_igblast7_vreg =~ s/#LC/$localconf/;   
    } else {
        $out_igblast7_vreg =~ s/#LC//;   
   }
if($out_igblast7_vreg eq $blast_result){
    print " !! ERROR: unexpected input file. .igblast7.tidy.txt expected !!\n";
    exit;
} else {
#    print "> peptide level report (best of blast7.tidy annotated) written to $out_igblast7_vreg.\n";
}

# 4. Other variables

my %ERRORS = ();
my @vregions = ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3");
my %vregion_template = ();
    foreach my $reg ((@vregions, "UNKNOWN")){ $vregion_template{$reg} = 0}


#  _              _ _
# | |___  __ _ __| (_)_ _  __ _
# | / _ \/ _` / _` | | ' \/ _` |
# |_\___/\__,_\__,_|_|_||_\__, |
#                         |___/
# 

# 1. INQ -> index accnr (=initial peptides sequence with version) and INQ peptide sequence
print "  > read accnr~pepseq index from $INQ_file ...";
my %INQ = ();
   %INQ = &read_INQ_index($INQ_file) if $useINQ == 1;
print scalar(keys(%INQ))." entries read.\n";


print "  > Reading regions:\n";
my %SEQ = &read_IMGT_fasta($IMGT_fasta_file);

# 2. Regions from IMGT IgBLAST

print "  > Loading IMGT V-gene regions from '$igblast7_tidy'..";
open(REG, $igblast7_tidy) || die "cannot open germline ig region files!/n";

my $reg_log = $igblast7_tidy;
   $reg_log =~ s/\.igblast7\.tidy\.txt/\.igblast7\.tidy\.log/;

if($reg_log eq $igblast7_tidy){
    print "!! UNexpected extension of v-region blast file. Ending with .igblast7.tidy.txt expected !!\n";
    exit;
}

my @REGLOG = ();

my $header  = <REG>;
chomp($header);
unless($header =~  /query.*alignment.*from.*to.*length.*matches.*mismatches.*gaps.*percent_identity.*subject_id/){
    print "!! Invalid IMGT-V-seq Igblastt tidy file !! (unexpected scan header)\n";
    print "   use parse-igblast7.pl on igblast file (format=7) to prepare file.\n";
    exit;
}

my %LEN = ();
my %POS = ();
my %REG = ();
my %GEN = ();
my %LAB = ();

my $rcnt = 0;
my $pcnt = 0;

while(<REG>){
    chomp;
    my ($query, $alignment, $from, $to, $length, $matches, $mismatches, $gaps, $percent_identity, $subject_id) = split(/\t/, $_);
    # query: M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | |98 AA|98+0=98| | |
    # see <https://www.imgt.org/IMGTindex/Fasta.php>
    my ($accnr, $gene, $tax, $allel_funct, $label, $range_na, $len_na, undef,undef,undef,undef, $len_aa, $n_char, undef ) = split(/\|/, $query);
    my $imgtid = join("|", ($accnr, $gene, $label));
    unless($n_char =~ /\+0=/){ print " !! unexpected term $n_char -> needs to be implemented !!\n"; exit;}

    if($alignment =~ /.*-IMGT/){
        ($alignment) = $alignment =~ /(.*)-IMGT/;
         $alignment =~ s/FR/FWR/;
        foreach my $pos ($from..$to){ $POS{$imgtid}{$pos} = $alignment; $pcnt++;}
        $REG{$imgtid}{$alignment}{start} = $from;
        $REG{$imgtid}{$alignment}{end}   = $to;

        if( defined $LEN{$imgtid}){
            $LEN{$imgtid} = $to if $LEN{$imgtid} < $to;
        } else {
            $LEN{$imgtid} = $to;
        }
        
        $rcnt++;
    }

    $GEN{$imgtid} = $gene;
    $LAB{$imgtid} = $label;
}

print ". done (".scalar(keys(%GEN))." ids, $rcnt regions and $pcnt postions).\n";

foreach my $imgtid (sort keys(%POS)){
    my %check_regions = %vregion_template;
    my %check_position = ();
    foreach my $p (1..$LEN{$imgtid}){ $check_position{$p} = 1}
    foreach my $p (sort keys(%{$POS{$imgtid}}) ) {
        my $r = $POS{$imgtid}{$p};
        if(defined $check_position{$p}){
            delete($check_position{$p});
        } else {
            print "!! $imgtid: more than one annotation on position $p !\n";
            exit;
        }
        if(defined $check_regions{$r}){delete($check_regions{$r})}
    }
    delete( $check_regions{"UNKNOWN"} );
    if(scalar(keys(%check_position)) > 0){
        #print "  > check $imgtid : ".scalar(keys(%check_position))." unannotated positions: ".join(" ", sort {$a <=> $b} keys(%check_position))."\n";
        push(@REGLOG, "  > check $imgtid : ".scalar(keys(%check_position))." unannotated positions: ".join(" ", sort {$a <=> $b} keys(%check_position)));
    }
    if(scalar(keys(%check_regions)) > 0){
        #print "  > check $imgtid : ".scalar(keys(%check_regions))." missing v-regions: ".join(" ", sort {$a cmp $b} keys(%check_regions))."\n";
        push(@REGLOG, "  > check $imgtid : ".scalar(keys(%check_regions))." missing v-regions: ".join(" ", sort {$a cmp $b} keys(%check_regions)));
    }
}

if(scalar(@REGLOG) > 0){
    print " > Writing ".scalar(@REGLOG)." warnings to $reg_log."; 
    open(REGLOG, ">$reg_log") || die $!;
    print REGLOG join("\n", @REGLOG)."\n";
    print "... done.\n";

}

close REG;

#                 _        _
#      _ __  __ _(_)_ _   | |___  ___ _ __
#     | '  \/ _` | | ' \  | / _ \/ _ \ '_ \
# ___ |_|_|_\__,_|_|_||_| |_\___/\___/ .__/
#                                    |_|


print "+ Start reading BLAST results ($igblast7_tidy).\n";
print "  - mode: ";
if($query_mode == 1){print "query mode (no output fileswill be  written).\n"} else { print "conversion mode (including writing of output files).\n"}
# sleep(1);
open(BLAST, $blast_result) || die "!! Cannot read BLAST result `$blast_result` !!\n";

my @OUT_PEPSEQ = ();
my @OUT_REGION = ();
my %OUT_BESTID = ();

my $proc_stp = 10000;
my $proc_cnt = 0;

print "> Running annotation ";

while(<BLAST>){

    chomp;
    # Read BLAST results (all seq are INQ sequences)

#print "\n\n   ($_)\n ";
    my @L = split(/\t/, $_);
    my ($qseqid, $sseqid, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $sseq, $qseq) =  split(/\t/, $_);
    my $imgtid = $sseqid;
    my ($accnr, $gene, $label) = split(/\|/, $sseqid);

    my $show_query_result = 0;

    if($query_mode == 1){
        if(scalar(@query_qseqid) > 0 && scalar(@query_imgtid) > 0){
            foreach my $q_qseqid (@query_qseqid){
                foreach my $q_imgtid (@query_imgtid){
                    # THIS NEEDS TO BE FIXED !!!!!
                    $q_imgtid =~ s/\*/x/g;
                    $q_imgtid =~ s/\|/ /g;
                    my $c_imgtid = $imgtid;
                       $c_imgtid =~ s/\*/x/g;
                       $c_imgtid =~ s/\|/ /g;
                    if($qseqid =~ /$q_qseqid/ && $c_imgtid =~ /$q_imgtid/){$show_query_result = 1}
                }
            }
        } elsif (scalar(@query_qseqid) > 0){
            foreach my $q_qseqid (@query_qseqid){
                    if($qseqid =~ /$q_qseqid/ ){ $show_query_result = 1 }
            }

        } elsif (scalar(@query_imgtid) > 0){
            foreach my $q_imgtid (@query_imgtid){
                    my $curr_imgtid = $imgtid;
                       $curr_imgtid =~ s/\*/x/g;
                       $curr_imgtid =~ s/\|/ /g;
                    if($curr_imgtid =~ /$q_imgtid/ ){ $show_query_result = 1 }
                    print "$imgtid (= $curr_imgtid) ~ $q_imgtid -> $show_query_result\n" if $show_query_result == 1; 
            }



        } else {
            print "! SOMEHTING WENT WRONG !\n\n";
            exit;
        }
    }


    if($query_mode == 0 || $show_query_result == 1){

            my ($pepseq);
            if($useINQ == 1){
                if(!defined $INQ{$qseqid}){
                    print "!! accession number $qseqid missing $INQ_file !!\n";
                } else {
                    $pepseq = $INQ{$qseqid};
                }
            } else {
                $pepseq = $qseqid;
            }

            my $qseq_align = $qseq;
               $qseq_align =~ s/-//g;

            my($unalign_n, $match, $unalign_c);
            my($len_unalign_n, $len_unalign_c);
            my($pstart, $pend);                    

            if($pepseq =~ $qseq_align){
                ($unalign_n, $match, $unalign_c) = $pepseq =~ /(.*)($qseq_align)(.*)/;
                $len_unalign_n = length($unalign_n);
                $len_unalign_c = length($unalign_c);
                $pstart        = $sstart - $len_unalign_n;
                $pend          = $send   + $len_unalign_c;
                my $sgaps    =  $sseq =~ tr/-//;
                my $qgaps    =  $qseq =~ tr/-//;
                my $calc_len = $pend - $pstart - $qgaps + $sgaps + 1;

                if(length($pepseq) != $calc_len ){
                    print "!! SOMETHING WHEN WRONG: peptide sequence != length from pstart and pend !!\n";
                    print " + qseq    : $qseq   -> $qseq --> length = ".length($qseq_align)."\n";
                    print "   * qgaps = $qgaps\n";
                    print " + sseq    : $sseq   -> length = ".length($sseq)."\n";
                    print "   * sgaps = $sgaps\n";
                    print " + sstart        : $sstart\n";
                    print " + send          : $send  \n";
                    print " + pepseq  : $pepseq -> length = ".length($pepseq)."\n";
                    print " + ".join(" | ",  ($unalign_n, $match, $unalign_c))."\n";
                    print " + len_unalign_n : $len_unalign_n\n";
                    print " + pstart        : $pstart\n";
                    print " + pend        : $pend\n";
                    print " + len_unalign_c : $len_unalign_c\n";
                    print " ! cal len : $calc_len\n\n";

                    print " + pend   : $pend\n";
                    exit if length($pepseq) == 12 && $qseq !~ /-/;
                } 


            } else {
                print " ! $pepseq NOT IN $qseq ($qseq_align) !\n";
                exit;
            }

            my $palign = sprintf("%.1f", length($match)*100/length($pepseq));

            if( $use_tags == 0 ){
              $TAG{$pepseq}{$pepseq} = -1;
            } elsif (!defined $TAG{$pepseq}){
                $TAG{$pepseq}{$qseq_align} = 0;
                print "!";
            }

            foreach my $t(keys(%{$TAG{$pepseq}}) ){

                my $tag_status = $TAG{$pepseq}{$t};
                my $tagseq     = $t;

                my ($lowconf_n, $conf, $lowconf_c) = $pepseq =~ /(.*)($t)(.*)/;
                my $len_lowconf_n = length($lowconf_n);
                my $len_lowconf_c = length($lowconf_c);

                # low/high confidence start/stops on query string
                my $qstartC = $len_lowconf_n + 1;
                my $qendC   = length($pepseq) - $len_lowconf_c;

                # low/high confidence start/stops on germine AA sequence
                my $pstartC = $pstart  + $len_lowconf_n;
                my $pendC   = $pend    - $len_lowconf_c;

                my $conf_pattern =  "."x$len_lowconf_n . ":"x(length($t)) . "."x$len_lowconf_c;
                my $pepseqCC_start = min($len_unalign_n, $len_lowconf_n);
                my $pepseqCC_end   = length($pepseq) - 1 - min($len_unalign_c, $len_lowconf_c);
                my $pepseqCC  = substr($pepseq, $pepseqCC_start, $pepseqCC_end - $pepseqCC_start + 1);
                my $pstartCC = $pstart + min($len_unalign_n, $len_lowconf_n);
                my $pendCC   = $pend   - min($len_unalign_c, $len_lowconf_c);

                my $report_id = join("@", ($pepseq, $label));


                my $report_flag = 0;

                if(!defined $OUT_BESTID{$report_id} || $evalue <= $OUT_BESTID{$report_id} ){
                    $OUT_BESTID{$report_id} = $evalue;
                    my $out_line = join("\t", ($qseqid, $sseqid, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $sseq, $qseq, $unalign_n, $unalign_c, $palign, $gene, $label, $pepseq, $tag_status, $tagseq, $pepseqCC, $qstartC, $qendC));
                    push(@OUT_PEPSEQ, $out_line);

                    $report_flag = 1;
                } 

                if($report_flag == 1){
                    if (($query_mode == 1 && $show_query_result == 1) || $debug == 1){
                        print  "\n\n__QUERY", $imgtid," & ", $qseqid,"_"x60, "\n";
                        print " + pepseq: ".$pepseq," aligned by $qseq_align [ => ", join(" | ", ($unalign_n, $match, $unalign_c))."]\n";
                        print " + tag sequence (local confidence > $localconf): $tagseq\n";
                        print " + alignment = ".$palign."% (aligned sequences / peptide sequence length)\n";
                    }

                    my @report_remark = ();

                    my $cstart   = $pstart;
                    my $cstartC  = $pstartCC;

                    if ($cstart < 1){
                        push(@report_remark, "coverage start <0 (".-$cstart.") -> set to 0");
                        $cstart = 1;
                    }
                    if ($cstartC < 1){
                        push(@report_remark, "coverage high conf start <0 (".-$cstart.") -> set to 0");
                        $cstartC = 1;
                    }

                    my $cend    = $pend;
                    my $cendC   = $pendCC;

                    if ($use_tags == 1){
                        $cstart = $pstartCC;
                        $cend   = $pendCC;
                    }

                    if(exists $POS{$imgtid}{$sstart}){ # thus, if == V-REGION !!!
                        my %curr_reg_a = %vregion_template;

                        foreach my $pos ($sstart..$send){
                            my ($pos_reg);
                            if(exists $POS{$imgtid}{$pos}){
                                $pos_reg = $POS{$imgtid}{$pos};
                            } else {
                                $pos_reg = "UNKNOWN";
                                printf ("[ %09d ]", $proc_cnt);
                                print " Unknown region for $imgtid at pos $pos -> set to UNKNOWN!\n";
                            }
                            $curr_reg_a{$pos_reg}++;
                        }

                        my %curr_reg_c = %vregion_template;


                        foreach my $pos ($cstart..$cend){
                            my ($pos_reg);
                            if(exists $POS{$imgtid}{$pos}){
                                $pos_reg = $POS{$imgtid}{$pos};
                            } elsif ($pos > $LEN{$imgtid}){
                                # all regions extending the n-terminus of the V-region are classified as CDR
                                $pos_reg = "CDR3";
                            } else {
                                $pos_reg = "UNKNOWN";
                                #printf ("[ %09d ]", $proc_cnt);
                                #print " Undefined region '$pos_reg' for $imgtid at pos $pos -> set to UNKNOWN!\n";
                                #print " not in: ".join(", ", keys(%vregion_template))."\n";
                            }
                            if(defined $curr_reg_c{$pos_reg}){
                                $curr_reg_c{$pos_reg}++;
                                
                            } else {
                                $curr_reg_c{$pos_reg} = 1;
                                print "! UNEXPECTED REGION $pos_reg !";
                                printf ("[ %09d ]\n", $proc_cnt);
                                print "  ~> $_\n";
                                print "  ~> current entries:\n";
                                foreach my $test_reg (sort keys(%curr_reg_c)){ print "\t -- $test_reg : ".$curr_reg_c{$test_reg }."\n";}
                                exit;
                            }
                        }

                        if($cend > $LEN{$imgtid}){
                            push(@report_remark, "c-term unaligned exceeds germline (",$cend - $LEN{$imgtid}," exceeding)");
                        }

                        if(($query_mode == 1 && $show_query_result == 1) || $debug == 1){
                            print "\n";
                            print "   - regions", " "x3;
                            foreach my $vreg (@vregions){
                                if(defined $REG{$imgtid}{$vreg}{end} && $REG{$imgtid}{$vreg}{start}){
                                    my $rlen = $REG{$imgtid}{$vreg}{end} - $REG{$imgtid}{$vreg}{start} + 1;
                                    my $plot_pre  = sprintf("%d", ($rlen - (2+length($vreg)))/2);
                                    my $plot_post = $rlen - (2+length($vreg)) - $plot_pre ;
                                    if($rlen >= 2+length($vreg)){
                                        print "<" , "~"x$plot_pre, $vreg, "~"x$plot_post, ">";
                                    } else {
                                        print "<" , "~"x($rlen - 2), ">";
                                    }
                                }
                            }

                            print "\n";
                            print "   - sequence", " "x2, $SEQ{$imgtid}." ($gene)\n";
                            print "   - sseq  ", " "x(4 + $sstart - 1), $sseq."\n";
                            print "   - match ", " "x(4 + $sstart - $len_unalign_n - 1 + $qstart - 1), &make_match_pattern($sseq, $qseq)."\n";
                            print "   - qseq  ", " "x(4 + $sstart - $len_unalign_n - 1 + $qstart - 1), $qseq."\n";
                            #print "   - pepseq", " "x(4 + $sstart - $len_unalign_n - 1), $pepseq."\n";
                            print "   - pepseq", " "x(4 + $pstart - 1), $pepseq."\n";
                            print "   - Locconf", " "x(3 + $pstart - 1), $conf_pattern."\n";
                            #pragseqint "   - tag   ", " "x(4 + $sstart - $len_unalign_n - 1 + $len_lowconf_n),$t."\n";
                            print "   - tag   ", " "x(4 + $pstartC - 1),$tagseq."\n";
                            print "   - confseq"," "x(3 + $pstartCC - 1), $pepseqCC."\n";
                            print "\n";

                        }
                        # ALGINED SEQUENCE
                        foreach my $vreg (@vregions){
                            if(exists $REG{$imgtid}{$vreg}{end} && exists $REG{$imgtid}{$vreg}{start}){
                                my $rlen = $REG{$imgtid}{$vreg}{end} - $REG{$imgtid}{$vreg}{start} + 1;
                                my $plen = length($pepseq);
                                my $slen = length($qseq);

                                my $reg_cov_seq = sprintf("%.1f", $curr_reg_a{$vreg}*100/$rlen);
                                my $seq_cov_reg = sprintf("%.1f", $curr_reg_a{$vreg}*100/$slen);

                                my $reg_cov_pep = sprintf("%.1f", $curr_reg_c{$vreg}*100/$rlen);
                                my $pep_cov_reg = sprintf("%.1f", $curr_reg_c{$vreg}*100/$plen);

                                if($reg_cov_pep > 0){
                                    # Query mode: 
                                    if (($query_mode == 1 && $show_query_result == 1) || $debug == 1){
                                        print "   * ALIGNED (seq) $vreg [".$REG{$imgtid}{$vreg}{start}."..".$REG{$imgtid}{$vreg}{end}."]: ".$curr_reg_a{$vreg}." aa -> $reg_cov_seq \%region covered \& $seq_cov_reg \%seq.\n";
                                        print "     COVERED (pep) $vreg [".$REG{$imgtid}{$vreg}{start}."..".$REG{$imgtid}{$vreg}{end}."]: ".$curr_reg_c{$vreg}." aa -> $reg_cov_pep \%region covered \& $pep_cov_reg \%seq.\n";
                                        print "\n";

                                    } 
                                    my ($out_remark);
                                    if(scalar(@report_remark) > 0){
                                        $out_remark = join("; ", @report_remark);
                                    } else {
                                        $out_remark = "NA";
                                    }
                                    my $out_line = join("\t", ($qseqid, $sseqid, $vreg, $curr_reg_c{$vreg}, $curr_reg_a{$vreg}, $rlen, $plen, $slen, $reg_cov_pep, $reg_cov_seq, $seq_cov_reg, $pep_cov_reg,  $out_remark, $pepseq, $tag_status, $t) ); # PEPSEQ AND TAG SEQ ARE INQSEQ!
                                    push(@OUT_REGION, $out_line);
                                } 
                            } else {
                                $ERRORS{$imgtid}++;

                            }
                        }
                        if($debug == 1){
                            print "DEBUG == $debug\n";
                            print "qseq         : $qseq\n";
                            print "sseq         : $sseq\n";
                            print "pepseq       : $pepseq\n";
                            print "tagseq       : $tagseq\n";
                            print "conf_pattern : $conf_pattern\n";
                            print "qstart       : $qstart\n";
                            print "qend         : $qend\n";
                            print "qstartC      : $qstartC\n";
                            print "qendC        : $qendC\n";
                            print "sstart       : $sstart\n";
                            print "send         : $send\n";
                            print " * len_unalign_n = $len_unalign_n\n";
                            print "pstart       : $pstart\n";
                            print " * len_unalign_c = $len_unalign_c\n";
                            print "pend         : $pend\n";
                            print "pstartC      : $pstartC\n";
                            print "pendC        : $pendC\n";
                            print "pstartCC     : $pstartCC\n";
                            print " * len_lowconf_n = $len_lowconf_n\n";
                            print "pendCC       : $pendCC\n";
                            print " * len_lowconf_c = $len_lowconf_c\n";

                            my $check_pepseqCC  = substr($pepseq, min($len_unalign_n, $len_lowconf_n), length($pepseq) -1 - min($len_unalign_c, $len_lowconf_c));
                            my $check_min_n     = min($len_unalign_n, $len_lowconf_n);
                            my $check_min_c     = min($len_unalign_c, $len_lowconf_c);
                            print " * check_min_n = $check_min_n (min of $len_unalign_n and $len_lowconf_n)\n";
                            print " * check_min_c = $check_min_c (min of $len_unalign_c and $len_lowconf_c)\n";
                            my $check_start = $check_min_n ;
                            print " * substr start  = check_min_n  = ".$check_start."\n";
                            my $check_end   = length($pepseq) - 1 - $check_min_c;
                            my $check_length  = $check_end - $check_start + 1;
                            print " * substr end    = length(".length($pepseq).")check_min_n - 1 - $check_min_n = ".$check_end."\n";
                            print " * substr length = end - start + 1 = ".$check_end."\n";
                            print " ~~> confident peptides : ".substr($pepseq, $check_start, $check_length)."\n";
                            <STDIN>;

                        }

                    } elsif ($label =~ /CH/ || $label =~ /C-REGION/ ){
                    # do nothing
                    } elsif($label =~ /J-REGION/ || $label =~ /D-REGION/) {

                        my $vreg = "CDR3".lc(substr($label,0,1));


                       if(($query_mode == 1 && $show_query_result == 1) || $debug == 1){
                            print "\n";
                            print "   - sequence", " "x22, $SEQ{$imgtid}." ($gene)\n";
                            print "   - sseq    ", " "x(22+$sstart-1), $sseq."\n";
                            print "   - pepseq "," "x(23 + $sstart - $len_unalign_n - 1), $pepseq."\n";
                            print "   - tag    "," "x(23 + $sstart - $len_unalign_n - 1 + $len_lowconf_n),$t."\n";
                            print "   - confseq"," "x(23 + $cstartC - 1), $pepseqCC."\n";
                            print "   - qseq   "," "x(23 + $sstart - $len_unalign_n - 1 + $qstart - 1), $qseq."\n";
                            print "   - match  "," "x(23 + $sstart - $len_unalign_n - 1 + $qstart - 1), &make_match_pattern($sseq, $qseq)."\n\n";
                            print "\n";
                        }
                       if(($query_mode == 1 && $show_query_result == 1) || $debug == 1){
                            print "\n";
                            print "   - sequence", " "x22, $SEQ{$imgtid}." ($gene)\n";
                            print "   - sseq  ", " "x(24 + $sstart - 1), $sseq."\n";
                            print "   - match ", " "x(24 + $sstart - $len_unalign_n - 1 + $qstart - 1), &make_match_pattern($sseq, $qseq)."\n";
                            print "   - qseq  ", " "x(24 + $sstart - $len_unalign_n - 1 + $qstart - 1), $qseq."\n";
                            print "   - pepseq", " "x(24 + $pstart - 1), $pepseq."\n";
                            print "   - Locconf", " "x(23 + $pstart - 1), $conf_pattern."\n";
                            print "   - tag   ", " "x(24 + $pstartC - 1),$tagseq."\n";
                            print "   - confseq"," "x(23 + $pstartCC - 1), $pepseqCC."\n";
                            print "\n";
                        }

                        my $rlen = length($SEQ{$imgtid});
                        my $plen = length($pepseq);
                        my $slen = length($qseq);
                        my $plenCC = length($pepseqCC);

                        my ($curr_reg_c);
                        if($use_tags == 1){
                            $curr_reg_c = length($pepseqCC);
                        } else {
                            $curr_reg_c = $slen;
                            $curr_reg_c += $len_unalign_n;
                        }

                        my $reg_cov_seq = sprintf("%.1f", $slen*100/$rlen);
                        my $reg_cov_pep = sprintf("%.1f", $curr_reg_c*100/$rlen);  # > 100% possible, because just part of the CDR3 region on D or J.
                        my $seq_cov_reg = sprintf("%.1f", $rlen*100/$plen);
                        my $pep_cov_reg = sprintf("%.1f", $curr_reg_c*100/$plen);
                        
                        if ($query_mode == 1 && $show_query_result == 1){
                            print "   * ALIGNED (seq) $vreg [".$label."]: ".$slen      ." aa -> $reg_cov_seq \%region covered \& $seq_cov_reg \%seq.\n";
                            print "     COVERED (pep) $vreg [".$label."]: ".$curr_reg_c." aa -> $reg_cov_pep \%region covered \& $pep_cov_reg \%seq.\n";
                            print "(usetags = $use_tags)\n";
                        } 

                        if($debug == 1){
                            print "DEBUG == $debug\n";
                            print "qseq         : $qseq\n";
                            print "sseq         : $sseq\n";
                            print "pepseq       : $pepseq\n";
                            print "tagseq       : $tagseq\n";
                            print "conf_pattern : $conf_pattern\n";
                            print "qstart       : $qstart\n";
                            print "qend         : $qend\n";
                            print "qstartC      : $qstartC\n";
                            print "qendC        : $qendC\n";
                            print "sstart       : $sstart\n";
                            print "send         : $send\n";
                            print " * len_unalign_n = $len_unalign_n\n";
                            print "pstart       : $pstart\n";
                            print " * len_unalign_c = $len_unalign_c\n";
                            print "pend         : $pend\n";
                            print "pstartC      : $pstartC\n";
                            print "pendC        : $pendC\n";
                            print "pstartCC     : $pstartCC\n";
                            print " * len_lowconf_n = $len_lowconf_n\n";
                            print "pendCC       : $pendCC\n";
                            print " * len_lowconf_c = $len_lowconf_c\n";

                            my $check_pepseqCC  = substr($pepseq, min($len_unalign_n, $len_lowconf_n), length($pepseq) -1 - min($len_unalign_c, $len_lowconf_c));
                            my $check_min_n     = min($len_unalign_n, $len_lowconf_n);
                            my $check_min_c     = min($len_unalign_c, $len_lowconf_c);
                            print " * check_min_n = $check_min_n (min of $len_unalign_n and $len_lowconf_n)\n";
                            print " * check_min_c = $check_min_c (min of $len_unalign_c and $len_lowconf_c)\n";
                            my $check_start = $check_min_n ;
                            print " * substr start  = check_min_n  = ".$check_start."\n";
                            my $check_end   = length($pepseq) - 1 - $check_min_c;
                            my $check_length  = $check_end - $check_start + 1;
                            print " * substr end    = length(".length($pepseq).")check_min_n - 1 - $check_min_n = ".$check_end."\n";
                            print " * substr length = end - start + 1 = ".$check_end."\n";
                            print " ~~> confident peptides : ".substr($pepseq, $check_start, $check_length)."\n";
                            <STDIN>;

                        }

                        my $out_line = join("\t", ($qseqid , $sseqid , $vreg , $curr_reg_c        , $slen              , $rlen , $plen , $slen , $reg_cov_pep , $reg_cov_seq , $seq_cov_reg , $pep_cov_reg , "by ".$label , $pepseq , $tag_status , $t ));
                        push(@OUT_REGION, $out_line);


                    }
                } 
            }
        }

    $proc_cnt++;
    print "." if($proc_cnt % $proc_stp == 0);
}
print "done.\n";



print "+ ERRORS &WARNINGS:\n";
print "  [!] NO START and/OR STOP values for ".scalar(keys(%ERRORS))." IMGT V germline sequences.\n";
print "  [!] Check ".scalar(@REGLOG)." errors and warnings about germline regions in $reg_log!\n" if scalar(@REGLOG) > 0;
print "\n\n";
if($query_mode == 0){
    print "< Writing to $out_igblast7_anno ...";
    open(OUT, ">$out_igblast7_anno") || die $!;

    #                     $qseqid,   $sseqid, $pident , $length , $mismatch , $gapopen , $qstart , $qend , $sstart , $send , $evalue , $bitscore , $sseq , $qseq , $unalign_n , $unalign_c , $palign , $gene , $label , $pepseq , $tag_status , $tagseq , $pepseqCC , $qstartC , $qendC));

    print OUT join("\t", ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq", "qseq", "unalign_n", "unalign_c", "palign", "gene", "label", "pepseq", "tag_status", "tagseq", "pepseqCC", "qstartC", "qendC"))."\n";
    foreach my $out_line (@OUT_PEPSEQ){
        print OUT $out_line."\n";
    }
    close OUT;
    print "done.\n";

    print "< Writing to $out_igblast7_vreg ...";
    open(OUT, ">$out_igblast7_vreg") || die $!;

    print OUT join("\t", ("qseqid", "sseqid", "vreg", "aa_reg_covered"   , "aa_reg_aligned"   , "rlen", "plen", "slen", "reg_cov_pep", "reg_cov_seq", "seq_cov_reg", "pep_cov_reg", "remark"     , "pepseq", "tag_status", "tag"))."\n";

    foreach my $out_line (@OUT_REGION){
        print OUT $out_line."\n";
    }
    close OUT;
    print "done.\n";
} else {
    print ".done (query mode -> no output written)."
}
close BLAST;

## subs ........................................................................


sub count_feature {
     my %R = ();
     foreach (@_){
          if (defined $R{$_}){
               $R{"$_"}++;
          }else{
               $R{"$_"} = 1;
          }
     }

     return %R;
}


sub read_IMGT_fasta {
    #V00558|IGKV1D-16*02|V-REGION
    my $fasta = shift(@_);
    $fasta =~ s/\\/\//g;
    my %S = ();
    my ($curr_id);
    open(FASTA, $fasta) || die "cannot read fasta $fasta!\n";
    my $lcnt = 0;
    my %IMGT_region_name = ();
    while(<FASTA>){
        $lcnt++;
        chomp;

        if($_ =~ /^>/){
            $_ =~ s/^>//;
            my @curr_accnr = split(/\|/, $_);
               $curr_id = join("|", ($curr_accnr[0], $curr_accnr[1], $curr_accnr[4]));
            if(defined $S{$curr_id}){
                print "!! ACCESSION NUMBER $curr_id NOT A UNIQUE KEY !!\n";
                exit;
            }

            if(defined $IMGT_region_name{$curr_accnr[4]}){
                push(@{$IMGT_region_name{$curr_accnr[4]}}, $curr_accnr[1]);
            } else {
                @{$IMGT_region_name{$curr_accnr[4]}} =  $curr_accnr[1];
            }

        } else {
            if(defined $S{$curr_id}){
                $S{$curr_id} .= $_;
            } else {
                $S{$curr_id} = $_;
            }
        }
    }

    foreach my $r (sort keys(%IMGT_region_name)) {
        my @iggenes = @{$IMGT_region_name{$r}};
        print "    -- $r : ".scalar(@{$IMGT_region_name{$r}})." entries: ".$iggenes[0]." ..".$iggenes[$#iggenes]."\n";
    }

    return %S;
}

sub read_INQ_index {
    my $index =shift(@_);
    open(IDX, $index) || die "cannot read index $index!\n";
    my %I = ();
    <IDX>;
    while(<IDX>){
        chomp;
        my($accnr, $pepseq, undef) = split(/\t/, $_);
        $I{$accnr} = $pepseq;
    }
    return %I;

}

sub make_match_pattern {
    if(length($_[0]) == length($_[1])){
        my $match = "";
        my @SSEQ = split(//, $_[0]);
        my @QSEQ = split(//, $_[1]);
        foreach my $pos (0..$#QSEQ){
            if($SSEQ[$pos] eq $QSEQ[$pos]){
                $match .= "+";
            } else {
                $match .= "-";
            }
        }
        return $match;
    } else {
        return "ERROR: unequal sequence length!\n";
    }
}


sub usage {
        print "Usage: blast-anno-igreg.pl de_novo_peptides_INQ.blast~IMGT-VDJC_v210614.txt (optional: query1,query2,query3)\n";
        print " Remarks: \n";
        print " + 'de_novo_peptides_INQ.blast~IMGT-VDJC_v210614.txt' files are required.\n";
        print " + the 'de_novo_peptides_INQ.txt' file needs to be in the same folder.\n";
        print " + queryies with | are recognized as IMGT-ID\n";
}


sub max {
     my @array = sort {$a <=> $b} @_;
     return $array[$#array];    
}

sub min {
     my @array = sort {$a <=> $b} @_;
     return $array[0];    
}

