#!/usr/bin/perl -w
use strict;
###########################################################################
##
## *  [2010] - [2016] Translational Genomics Research Institute (TGen)
## *  All Rights Reserved.
## *
## * Major Contributor(s): David W. Craig, Ph.D.
##
##  This script produces a read count at every 100 bases across entire genome.
##  It leverages read pairs by considering the distance between read pairs as
##  being covered (physical coverage). If the position is covered by an insert(s)
##  (=read1+spanning distance+read2), it will have a read count associated with it.
##  Otherwise, the position will have a value of 0. Output file is a three column
##  tab-delimited text file with .dat extension. Column 1 is chromosome, column 2
##  is position, column 3 is read count.
##  Currently, this script is only implemented for the human genome. 
##  Chromesomes X,Y and MT are converted to 23, 24, and 25, respectively.
##
##  Input:
##	I= BAM from whole genome, exome or custome panel sequencing.
##  O= Filename for output DAT file. Extension '.dat' is automatically appended.
##	M= RG:  Optional matching of reads containing RG:
##	S= samtools path (including path to executable)
##
##  Output:
##	DAT file with three columns: Chromosome, Position and Read Count
##
##  Example:
##	./tgen_CloneCov.pl \
#		I=${BAMFILE} \
#	        O=${OUTFILE} \
#	        M=RG: \
#	        S=${SAMTOOLSPATH}/samtools
##
##  $Id: tgen_CloneCov.pl,v 1.01 dcraig Exp $
##
############################################################################

#######################
# Defaults
#######################
my $VERS  = 1.01;
my $bin   = 100;
my $flags = "-f 0x0002 -F 0x0400 -F0x0200 -q 20";

my @ln_chr = ();
my @cln    = ();
my $out   = "MISSING_Output_header";
my $ma    = "RG:";
my $st    = "samtools";
my $i_min=10;
my $i_max=2000;
my $mqf=20;
my $in = "MISSINGFILE.bam";
#######################
# Command Line Processing
#######################
print "Script:\ttgen_CloneCov.pl (version:$VERS):\tEstimates physical coverage to nearest 100bp\n";
print "\t-Usage:\t\t./tgen_clone_coverage.pl I=myInput.bam S=/my/samtoolsPath/samtools M=RG: O=outputName > outputLog\n";
print "\n\t-Arguments:\t@ARGV\n";
if ( $#ARGV < 1 || join(@ARGV) =~ /-h/ ) {
    print "\t-Not enough arguments, exiting\n";
    exit;
} else {
    foreach my $a (@ARGV) {
        if ( $a =~ /I=/ ) {
            $in = $a;
            $in =~ s/I=//g;
            print "\t-Input_file:\t$in\n";
        } elsif ( $a =~ /O=/ ) {
            $out = "$a.cln";
            $out =~ s/O=//g;
            print "\t-Output file:\t$out.cln.dat\n";
        }elsif ( $a =~ /S=/ ) {
            $st = $a;
            $st =~ s/S=//g;
            print "\t-Samtools path: $st\n";
        }elsif ( $a =~ /M=/ ) {
            $ma = $a;
            $ma =~ s/M=//g;
            print "\t-Match BAM records with: $ma\n";
        }else {
            die "DIE: Can't parse $a\n";
        }
    }
}
#######################
# Initialize Coverage
#######################
print "\n+Running and Initializing\n";
open( IN, "$st view -H $in |" ) || die "Can't open $in\n";
while (<IN>) {
    chomp;
    my @buf = split(/\t/);
    if ( $buf[1] =~ /SN:/ ) {
        my $chr = $buf[1];
        $chr =~ s/SN://g;
        $chr =~ s/chr//g;
        my $len = $buf[2];
        $len =~ s/LN://g;
        $len = int( $len / $bin ) + 1;
        if ( $chr eq "X" )  { $chr = 23 }
        if ( $chr eq "Y" )  { $chr = 24 }
        if ( $chr eq "MT" ) { $chr = 25 }

        if ( $chr>0 && $chr<=25 ) {
            $ln_chr[$chr] = $len;
            print "\t-Initializing chr: $chr to $len\n";
            for ( my $i = 0 ; $i <= $len ; ++$i ) {
                $cln[$chr][$i] = 0;
            }
        }
    }
}
close(IN);
print "+Done Initializing\n";
#################
# Main Loop
#################
my $c = 0;
open( IN, "$st view $flags $in |" ) || die "Can't open $in\n";
print "+Start Processing BAM $in\n";
LOOP: while (<IN>) {
    my @temp = split(/\t/);
    ++$c;
    if ( $c % 10000000 == 0 ) { print "\t-Read: $c at $temp[2]\t$temp[3]\n"; }
    $temp[2] =~ s/chr//g;
    if ( $temp[2] eq "X" )  { $temp[2] = 23; }
    if ( $temp[2] eq "Y" )  { $temp[2] = 24; }
    if ( $temp[2] eq "MT" ) { $temp[2] = 25; }
    if ( $temp[2] =~ /\D/ ) { 
       next LOOP; 
    }
	my $chr = $temp[2];    
    if ( $chr>0 && $chr<=25 ) {
		if ( $temp[8] > $i_min && $temp[8] < $i_max ) {
			my $st  = int( $temp[3] / $bin );
			my $en  = int( $temp[7] / $bin );
			for ( my $j = $st ; $j <= $en ; ++$j ) {
				++$cln[$chr][$j];
			}
			next LOOP;
		}
    }
}
close(IN);
print "\t-Reads Processed: $c\n";
#################
# Print Output
#################
print "+Printing output\n";
open( OUT, ">$out.dat" ) || die "Can't open $out.dat\n";
for ( my $i = 1 ; $i <= 25 ; ++$i ) {
    for ( my $j = 0 ; $j <= $ln_chr[$i] ; ++$j ) {
        my $pos = $j * $bin;
		if ( exists( $cln[$i][$j] ) ) {
			print OUT "$i\t$pos\t$cln[$i][$j]\n";
		}
    }
}
close(OUT);
print "Done.\n";
