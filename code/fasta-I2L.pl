#!/usr/bin/perl -w

use strict;
$| = 1;

#print "\n" x 2;
print "+ replace Ile -> Leu (Christoph Stingl, 17. May 2022)\n\n";

foreach ( 0 .. $#ARGV ) { $ARGV[$_] =~ s/\\/\//g; }

if(!defined $ARGV[0]){
    print "usage: fasta-I2L.pl input1.fasta";
    exit;
} 

my $fasta_file = shift(@ARGV);
my $ILout_file = $fasta_file;
   $ILout_file =~ s/\.fasta/\_I2L\.fasta/i;

if($fasta_file eq $ILout_file){
    print "!! No fasta file specified !!\n";
    print "usage: fasta-to-tbl.pl input1.fasta";

    exit;
}

my %AA = ();

unless ( -e $fasta_file ) {
    print "!! $fasta_file does not exists!!\n";
    print "usage: fasta-to-tbl.pl input1.fasta";

    exit;
} else {

    print "< reading $fasta_file\n";

    my %entry_lines   = ();
    my $entry_line    = 0;
    my $current_start = 0;
    my $cnt_e = 0;
    my $cnt_l = 0;
    my $cnt_I = 0;
    my $cnt_aa = 0;
    open( IN, $fasta_file ) || die "Can not read sfd file 'fasta_input_file'.\n";
    open( IL, ">$ILout_file" ) || die "Can not write to new fasta file.\n";

    while (<IN>) {
        chomp;

        $cnt_l++;
        if ( $_ =~ /^>/ ) {
            $cnt_e++;
        } else {
            my @S = split(//, $_);
            my %A = &count_feature(@S);
            foreach my $aa  (keys (%A)){
                if(defined $AA{$aa}){
                    $AA{$aa} += $A{$aa};
                } else {
                    $AA{$aa} = 1;
                }
            }
            $_ =~ s/I/L/g;
            $cnt_aa += length($_);
        }
        print IL $_."\n";
    }

    close IN;
    close IL;

    print "SUMMARY";
    print "_"x73,"\n";
    print "1.) entries        : $cnt_e\n";
    print "2.) lines          : $cnt_l\n";
    print "3.) aminoacids     : $cnt_aa\n";
    print "4.) All amino acids:\n";
    foreach my $aa (sort keys(%AA)){
        print "    + n $aa: ".$AA{$aa}."\n";
    }

}


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