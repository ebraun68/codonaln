#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use POSIX;

############################################################################
# checkframe.pl
#
#  - Check the frame for a set of nucleotide sequences and output a fasta
#    file with sequences padded with N's
############################################################################

my($progname) = $0;

############################################################################
# Initialize variables
############################################################################

my($iter);
my($jter);
my($kter);
my($lter);my($mter);
my($zter);

my($tval1);
my($tval2);
my($tstring);
my @temparray;

if ( @ARGV < 3 ) {
	print "Usage:\n  \$ $progname <infile> <outfile> <STOP1> <STOP2> <STOP3> <STOP4>\n";
	print "   The first two arguments are mandatory\n";
	print "      infile  = tab-delimited table of sequence names and sequences\n";
	print "      outfile = prefix for output files\n";
	print "      - This program will write two output files:\n";
	print "          <outfile>.bestframe.fasta  -- nucleotide sequences in correct frame\n";
	print "          <outfile>.frames.txt       -- numbers of stop codons in each frame\n";
	print "      STOPx   = list of stop codons to consider\n";
	print "      - Use STANDARD or STD for the standard genetic code (TAA, TAG, TGA)\n";
	print "      - Use VMT for the vertebrate mitochondrial code (TAA, TAG, AGA, AGG)\n";
	print "      - Otherwise list the stop codons to consider\n";
	print "          For example:\n";
	print "            TAA, TAG, TGA -- Standard code without ambiguities\n";
	print "            TGA           -- Ciliate, Dasycladacean, and Hexamita nuclear code\n";
	print "          Up to eight stop codons (including ambiguities) can be specified\n";
	print "exiting...\n";
	exit;
}

my($infile)  = $ARGV[0];
my($outfile) = $ARGV[1];
my($nstops)  = @ARGV - 2;
my @STOP;
if ( uc($ARGV[2]) eq "STANDARD" || uc($ARGV[2]) eq "STD" ) {
	$STOP[0] = "TAA";
	$STOP[1] = "TAG";
	$STOP[2] = "TGA";
	$STOP[3] = "TAR";
	$nstops = 4;
}
elsif ( uc($ARGV[2]) eq "VMT" ) { 
	$STOP[0] = "TAA";
	$STOP[1] = "TAG";
	$STOP[2] = "AGA";
	$STOP[3] = "AGG";
	$STOP[4] = "TAR";
	$STOP[5] = "AGR";
	$nstops = 6;
}
else { # alt code, stop codons passed on command line
	$STOP[0] = uc($ARGV[2]);
	if ( @ARGV > 3 ) { $STOP[1] = uc($ARGV[3]); }
	if ( @ARGV > 4 ) { $STOP[2] = uc($ARGV[4]); }
	if ( @ARGV > 5 ) { $STOP[3] = uc($ARGV[5]); }
	if ( @ARGV > 6 ) { $STOP[4] = uc($ARGV[6]); }
	if ( @ARGV > 7 ) { $STOP[5] = uc($ARGV[7]); }
	if ( @ARGV > 8 ) { $STOP[6] = uc($ARGV[8]); }
	if ( @ARGV > 9 ) { $STOP[7] = uc($ARGV[89]); }
}

print "Check sequence frames\n\n";
print "Consider for a total of $nstops stop codons:\n";
print "@STOP\n\n";

############################################################################
# Read the input file
############################################################################
print "Reading input file $infile... ";
open (my $INF, $infile) or die "Could not open file $infile for input.\n";
my @seqdata = <$INF>; # Read the input file
close($INF) or die "Could not close file $infile.\n";
my($taxa) = $#seqdata + 1;

my @name;
my @seq;

for ($iter=0; $iter<$taxa; $iter++) {
	($name[$iter],$seq[$iter]) = split(/\s+/, $seqdata[$iter]);
}

print "done.\n\nInput file has $taxa taxa\n\n";

############################################################################
# Check the number of stop codons in each frame
############################################################################

# Forward frame 1 = 0 -- Reverse frame 1 = 3
# Forward frame 2 = 1 -- Reverse frame 2 = 4
# Forward frame 3 = 2 -- Reverse frame 3 = 5
my @correctedseq;
my @inframestops;
my @incomplete_codons;

my($seqlen);
my($codon);

my($bestframe);	my($minstops);

my @frameID = ( "F1", "F2", "F3", "R1", "R2", "R3" );

open (my $DAF, ">$outfile.frames.txt") or die "Could not open file $outfile.frames.txt for output.\n";
open (my $FAF, ">$outfile.bestframe.fasta") or die "Could not open file $outfile.bestframe.fasta for output.\n";

print $DAF "Name\tF1stop\tF1gap\tF2stop\tF2gap\tF3stop\tF3gap\tR1stop\tR1gap\tR2stop\tR2gap\tR3stop\tR3gap\tBest_frame\tMinStops\n";

for ($iter=0; $iter<$taxa; $iter++) {
	print "$name[$iter]";
	print $DAF "$name[$iter]";
	print $FAF ">$name[$iter]";
	
	$correctedseq[0] = uc($seq[$iter]);
	$correctedseq[1] = "NN" . "$seq[$iter]";
	$correctedseq[2] = "N" . "$seq[$iter]";
	
	# correctedseq[3] is reverse complement frame 1
	$correctedseq[3] = reverse $seq[$iter];
	$correctedseq[3] =~ tr/ATGCatgc/TACGtacg/;
	# ambiguity codes written separately to improve readability of perl code
	$correctedseq[3] =~ tr/RYKMrykm/YRMKyrmk/;
	$correctedseq[3] =~ tr/BVDHbvdh/VBHDbvhd/;
	$correctedseq[4] = "NN" . "$correctedseq[3]";
	$correctedseq[5] = "N" . "$correctedseq[3]";
	
	# pad the 3' end of the sequences
	for ($jter=0; $jter<6; $jter++) {
		if ( length($correctedseq[$jter]) % 3 == 1 ) { $correctedseq[$jter] = "$correctedseq[$jter]" . "NN"; }
		if ( length($correctedseq[$jter]) % 3 == 2 ) { $correctedseq[$jter] = "$correctedseq[$jter]" . "N"; }
#		print "$jter - $correctedseq[$jter]\n\n";
	}
	
	# iterate over all 6 frames
	for ($jter=0; $jter<6; $jter++) {
	
		$seqlen = length($correctedseq[$jter]) / 3.0;
#		print "$seqlen - ";
		$inframestops[$jter] = 0;
		$incomplete_codons[$jter] = 0;
		for ($kter=0; $kter<$seqlen; $kter++) {
			$zter = $kter * 3;
			$codon = substr($correctedseq[$jter],$zter,3);
#			print "$codon ";
			# check for stop codons
			for ($lter=0; $lter<$nstops; $lter++) {
				if ( $codon eq $STOP[$lter] ) { $inframestops[$jter]++; }
				if ( $codon =~ /-/ && $codon ne "---" ) { $incomplete_codons[$jter]++; }
			}
		}
	
		if ( $jter == 0 ) { $bestframe=0; $minstops = $inframestops[$jter]; }
		elsif ( $inframestops[$jter] < $minstops ) { $bestframe = $jter; $minstops = $inframestops[$jter] }
		print $DAF "\t$inframestops[$jter]\t$incomplete_codons[$jter]";
	
	}
	
#	print "\t$bestframe";
	print "\t$frameID[$bestframe]\n";
	print $DAF "\t$frameID[$bestframe]\t$inframestops[$bestframe]\n";
	print $FAF "_$frameID[$bestframe]";
	print $FAF "_frame\n";
	print $FAF "$correctedseq[$bestframe]\n";
	
}

close($DAF) or die "Could not close file $outfile.frames.txt\n";
close($FAF) or die "Could not close file $outfile.bestframe.fasta\n";


