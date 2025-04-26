#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use POSIX;

############################################################################
# codon_align.pl
#
#  - Translate an unaligned fasta file of coding sequences to amino acids,
#    align the protein sequences, and use the amino acid alignment to guide
#    the nucleotide alignment (i.e., the codon alignment)
############################################################################

my($progname) = $0;

# paths to aligners (may need to change)
my($muscleexec)  = "muscle";
my($muscle5exec) = "muscle";
my($mafftexec)   = "mafft";

############################################################################
# Initialize variables
############################################################################

my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($zter);

my($tval1);
my($tval2);
my($tstring);
my @temparray;

if ( @ARGV < 2 || @ARGV > 6 ) {
	print "Usage:\n  \$ $progname <infile> <outfile> OPTIONAL: [aligner] [--protaln filename] [-M or --multiline]\n";
	print "   The first two arguments are mandatory\n";
	print "      infile  = single line aligned fasta infile\n";
	print "      outfile = prefix for output files\n";
	print "      - This program will write five output files:\n";
	print "          <outfile>.unaligned.faa         -- amino acid sequences\n";
	print "          <outfile>.aligned.faa           -- aligned amino acid sequences\n";
	print "          <outfile>.seq_summary.txt       -- summary of sequence lengths and errors\n";
	print "          <outfile>.codon.alignment.fasta -- codon aligned sequences (fasta format)\n";
	print "          <outfile>.codon.alignment.phy   -- codon aligned sequences (phylip format)\n";
	print "   The third and fourth arguments are optional\n";
	print "      - muscle (v3) is the default aligner (can also specified using --muscle3)\n";
	print "      - mafft can be specified using --mafft\n";
	print "      - muscle v5 align mode can be specified using --muscle5\n";
	print "      - muscle v5 super5 mode can be specified using --super5\n";
	print "      - pass --protaln followed by a filename to use an existing protein alignment \n";
	print "\n      - use --translate to translate sequences without aligning\n";
	print "      - use --check to assess numbers of in frame stop codons without aligning sequences\n";
	print "\n      - If -M or --multiline is passed the input file can be a multiline fasta file\n";
	print "      - The multiline fasta file will be converted to singleline and saved with the\n";
	print "        the filename infile.singleline\n";
	print "      - If you want to use a multiline input file but do not want to save a the singleline\n";
	print "        file use the clean multiline mode by passing -MC or --multiclean\n";
	print "exiting...\n";
	exit;
}

my($infile)  = $ARGV[0];
my($outfile) = $ARGV[1];
my($aligner) = "muscle";
my($multiline) = 0;
my($clean)     = 0;
my($protalnfile);
if ( @ARGV > 2 ) {
	if ( $ARGV[2] eq "--mafft" ) { $aligner="mafft"; }
	if ( $ARGV[2] eq "--muscle5" ) { $aligner="muscle5"; }
	if ( $ARGV[2] eq "--super5" ) { $aligner="super5"; }
	if ( $ARGV[2] eq "--protaln" ) { 
		$aligner="protaln";
		if ( @ARGV < 4 ) {
			print "ERROR: a protein alignment file name is necessary to run in --protaln mode\n";
			print "exiting...\n";
			exit;
		}
		$protalnfile = $ARGV[3];
		if ( ! -e $protalnfile ) {
			print "ERROR: protein alignment file $protalnfile does not exist\n";
			print "exiting...\n";
			exit;
		}
	}
	if ( $ARGV[2] eq "--check" ) { $aligner="check"; }
	if ( $ARGV[2] eq "--translate" ) { $aligner="translate"; }
	if ( $ARGV[2] eq "--multiline" || $ARGV[2] eq "-M" ) { $multiline=1; }
	if ( $ARGV[2] eq "--multiclean" || $ARGV[2] eq "-MC" ) { $multiline=1; $clean=1; }
}
if ( @ARGV > 3 ) {
	if ( $ARGV[3] eq "--mafft" ) { $aligner="mafft"; }
	if ( $ARGV[3] eq "--muscle5" ) { $aligner="muscle5"; }
	if ( $ARGV[3] eq "--super5" ) { $aligner="super5"; }
	if ( $ARGV[3] eq "--protaln" ) { 
		$aligner="protaln";
		if ( @ARGV < 5 ) {
			print "ERROR: a protein alignment file name is necessary to run in --protaln mode\n";
			print "exiting...\n";
			exit;
		}
		$protalnfile = $ARGV[4];
		if ( ! -e $protalnfile ) {
			print "ERROR: protein alignment file $protalnfile does not exist\n";
			print "exiting...\n";
			exit;
		}
	}
	if ( $ARGV[3] eq "--check" ) { $aligner="check"; }
	if ( $ARGV[3] eq "--translate" ) { $aligner="translate"; }
	if ( $ARGV[3] eq "--multiline" || $ARGV[3] eq "-M" ) { $multiline=1; }
	if ( $ARGV[3] eq "--multiclean" || $ARGV[3] eq "-MC" ) { $multiline=1; $clean=1; }
}
if ( @ARGV > 4 ) {
	# $ARGV[4] should either be the filename passed with --protaln (in which case the
	# variable is read above) or -M, --multiline, -MC, or --multiclean.
	if ( $ARGV[4] eq "--multiline" || $ARGV[3] eq "-M" ) { $multiline=1; }
	if ( $ARGV[4] eq "--multiclean" || $ARGV[3] eq "-MC" ) { $multiline=1; $clean=1; }
}

if ( $aligner eq "protaln" && "$protalnfile" eq "$outfile.aligned.faa" ) {
	print "ERROR: The alignment used in --protaln mode have the name of an output file\n";
	print "  -- in this case, the outfile name is $outfile.aligned.faa\n";
	print "exiting...\n";
	exit;
}

# Note that it is possible to pass multiple aligner names; this would be an error that
# the program does not check, although the run will tell the user which aligner has been
# chosen.


############################################################################
# Hash for genetic code (standard code)
#   NOTE: most third position ambiguities will yield the correct amino acid
#   NOTE: threefold ambiguities in third positions are ignored - will yield a "non-codon"
#   NOTE: stop codons set to X
# any non-codon will yield X because this construct is used:
#     $aa = exists($stdgencode{$codon}) ? $stdgencode{$codon} : 'X';
#
my %stdgencode = (
# Unambiguous codons
TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L', 
CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L', 
ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M', 
GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V', 

TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S', 
CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P', 
ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T', 
GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A', 

TAT => 'Y', TAC => 'Y', TAA => 'X', TAG => 'X', 
CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q', 
AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K', 
GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E', 

TGT => 'C', TGC => 'C', TGA => 'X', TGG => 'W', 
CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R', 
AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R', 
GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G', 

# Codons with ambiguities
TTY => 'F', TTR => 'L', 
CTY => 'L', CTR => 'L', CTS => 'L', CTW => 'L', CTK => 'L', CTM => 'L', CTN => 'L', YTR => 'L',
ATY => 'I', ATW => 'I', 
GTY => 'V', GTR => 'V', GTS => 'V', GTW => 'V', GTK => 'V', GTM => 'V', GTN => 'V', 

TCY => 'S', TCR => 'S', TCS => 'S', TCW => 'S', TCK => 'S', TCM => 'S', TCN => 'S', 
CCY => 'P', CCR => 'P', CCS => 'P', CCW => 'P', CCK => 'P', CCM => 'P', CCN => 'P', 
ACY => 'T', ACR => 'T', ACS => 'T', ACW => 'T', ACK => 'T', ACM => 'T', ACN => 'T', 
GCY => 'A', GCR => 'A', GCS => 'A', GCW => 'A', GCK => 'A', GCM => 'A', GCN => 'A', 

TAY => 'Y', 
CAY => 'H', CAR => 'Q', 
AAY => 'N', AAR => 'K', 
GAY => 'D', GAR => 'E', 

TGY => 'C', 
CGY => 'R', CGR => 'R', CGS => 'R', CGW => 'R', CGK => 'R', CGM => 'R', CGN => 'R', MGA => 'R', MGG => 'R', MGR => 'R',
AGY => 'S', AGR => 'R', 
GGY => 'G', GGR => 'G', GGS => 'G', GGW => 'G', GGK => 'G', GGM => 'G', GGN => 'G'
);

############################################################################
# Convert to single line fasta if --multiline flag is set
############################################################################
if ( $multiline == 1 ) {
	# Use awk to convert to singline. General version of the awk command is:
	# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < multiline.fasta | tail -n +2 > singleline.fasta
	system("awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"%s\",\$0);}  END {printf(\"\\n\");}' < $infile | tail -n +2 > $infile.singleline");
	$infile = "$infile" . ".singleline";
}

############################################################################
# Read the input file
############################################################################
print "Reading input file $infile... ";
open (my $INF, $infile) or die "Could not open file $infile for input.\n";
my @seqdata = <$INF>; # Read the input file
close($INF) or die "Could not close file $infile.\n";
my($seqlines) = $#seqdata + 1;

my($taxa) = 0;
my($maxnamelength) = 10;

for ($iter=0; $iter<$seqlines; $iter++) {
	if ( $seqdata[$iter] =~ /\>/ ) { 
		$taxa++;
		if ( length($seqdata[$iter]) > $maxnamelength ) { $maxnamelength = length($seqdata[$iter]); }
	}
}

print "done.\n\nInput file has $taxa taxa\n\n";

############################################################################
# Output a translated unaligned fasta file
############################################################################
my($aa);
my($len);
my($codon);

my @name;
my @seq;
my @numX;
my @sites;
my @frame;
my @stops;

print "Generating translated fasta file... ";

open (my $OUTF, ">$outfile.unaligned.faa") or die "Could not open file $outfile.unaligned.faa for output.\n";

for ($iter=0; $iter<$taxa; $iter++) {
	$jter = $iter * 2;
	$name[$iter] = $seqdata[$jter];
	chomp($name[$iter]);
	$name[$iter] =~ s/>//g;
	print $OUTF ">$name[$iter]\n";
	$seq[$iter] = $seqdata[$jter+1];
	chomp($seq[$iter]);
	$sites[$iter] = length($seq[$iter]);
	$len = $sites[$iter] / 3;
	$frame[$iter] = $sites[$iter] % 3;
#	will echo the length of sequences - commented out to reduce rubbish printed to screen
#	print "Length of sequence $name[$iter] is $sites[$iter] (modulo 3 = $frame[$iter])\n";
	# pad the sequence with N's if $frame !=0
	if ( $frame[$iter] == 1 ) { $seq[$iter] = "$seq[$iter]" . "NN"; }
	if ( $frame[$iter] == 2 ) { $seq[$iter] = "$seq[$iter]" . "N"; }
	# zero out the relevant element of @numX
	$numX[$iter] = 0;
	$stops[$iter] = 0;
	for ($kter=0; $kter<$len; $kter++) { 
		$mter = $kter * 3;
		$codon = substr($seq[$iter],$mter,3);
		$codon = uc($codon);
		$aa = exists($stdgencode{$codon}) ? $stdgencode{$codon} : 'X';
		if ( $aa eq "X" ) { $numX[$iter]++; }
		print $OUTF "$aa";
		# check for stop codons given the universal code
		# (change this if you desire another code)
		if ( $codon eq "TAA" || $codon eq "TAG" ) { $stops[$iter]++; }
		elsif ( $codon eq "TGA" ) { $stops[$iter]++; }
	}
#	print "   -- undefined amino acids = $numX[$iter]\n";
	print $OUTF "\n"; 
}

close($OUTF) or die "Could not close file $outfile.unaligned.faa\n";

print "done\n";

############################################################################
# Write tab-delimited file with sequence quality information. This file has
# four columns:
#    -- sequence name
#    -- number of nucleotides in the sequence (unaligned)
#    -- ending frame of the sequence (# of nucleotides modulo 3)
#    -- number of in frame stop codons
#    -- number of X's (undefined amino acids) in the translated sequence
# NOTE: this program assumes the first nucleotide of all sequences is the
#       first position of a codon
############################################################################
my($maxX)  = 0;
my($maxSC) = 0; # number of stop codons
my($numFP) = 0; # number of "frame problems"

print "\nPrinting summary of translations\n";
open (my $TABF, ">$outfile.seq_summary.txt") or die "Could not open file $outfile.seq_summary.txt for output.\n";

print $TABF "Name\tSites\tEndFrame\tNstops\tNumX\n";
for ($iter=0; $iter<$taxa; $iter++) {
	print $TABF "$name[$iter]\t$sites[$iter]\t$frame[$iter]\t$stops[$iter]\t$numX[$iter]\n";
	if ( $numX[$iter] > $maxX ) { $maxX = $numX[$iter]; }
	if ( $frame[$iter] > 0 ) { $numFP++; }
	if ( $stops[$iter] > $maxSC ) { $maxSC = $stops[$iter]; }
}

close($TABF) or die "Could not close file $outfile.seq_summary.txt\n";

print "\nTotal number of sequences                               = $taxa\n";
print "Number of sequences that are not a multiple of three    = $numFP\n";
print "Maximum number of in frame stop codons in any sequence  = $maxX\n";
print "Maximum number of undefined amino acids in any sequence = $maxX\n";

print "\nWARNING: Check sequences with an excessive number of in frame stop codons or undefined amino acids\n";

if ( $aligner eq "check" || $aligner eq "translate" ) {

	# Clean up the singleline fasta file if mode is --multiclean or -MC
	if ( $multiline == 1 && $clean == 1 ) { system("rm -f $infile"); }

	# Remove the unaligned output file an echo appropriate message if mode is --check
	if ( $aligner eq "check" ) { 
		system("rm -f $outfile.unaligned.faa"); 
		
		print "\nProgram was run in --check mode to assess sequence quality. Sequences will not be aligned.\n\n";
		print "NOTE: Running in --check mode only produces one output file:\n";
		print "         $outfile.seq_summary.txt\n";
		print "      This file can be used to assess the quality of each sequence\n";		
	}
	else { # echo appropriate message if mode is --translate
		print "\nProgram was run in --translate mode to generate protein sequence. Sequences will not be aligned.\n\n";
		print "NOTE: Running in --translate mode only produces one output file:\n";
		print "         $outfile.unaligned.faa\n";
		print "      This file can be used for other analyses (e.g., input for other aligners)\n";
	}

	print "\nexiting...\n\n";
	exit;
}

############################################################################
# Perform the sequence alignment
############################################################################
if ( $aligner eq "protaln" ) { print "\nCodon alignment will be generated using $protalnfile\n"; }
else { print "\nSequences will be aligned using $aligner\n\n"; }

# Remove the aligned output file if it exists
system("rm -f $outfile.aligned.faa");

# Additional alignment programs can be added here, but their ultimate output should be
# aligned fasta for this program to work
if ( $aligner eq "mafft" ) {
	system("$mafftexec --auto $outfile.unaligned.faa > $outfile.aligned.faa");
	
	# Use the following if absence of the aligner will generate an empty file
	if ( -s "$outfile.aligned.faa" == 0 ) {
		print "Alignment was not generated. Check path to mafft aligner\n\n";
		print "exiting...\n";
		system("rm -f $outfile.aligned.faa");
		exit;
	}
}
elsif ( $aligner eq "muscle5"  ) { 
	system("$muscle5exec -align $outfile.unaligned.faa -output $outfile.aligned.faa");
	
	unless ( -e "$outfile.aligned.faa" ) {
		print "Alignment was not generated. Check path to muscle aligner\n\n";
		print "exiting...\n";
		exit;
	}
}
elsif ( $aligner eq "super5" ) { 
	system("$muscle5exec -super5 $outfile.unaligned.faa -output $outfile.aligned.faa");
	
	unless ( -e "$outfile.aligned.faa" ) {
		print "Alignment was not generated. Check path to muscle aligner\n\n";
		print "exiting...\n";
		exit;
	}
}
elsif ( $aligner eq "protaln" ) {
	# copy the protein alignment file to $outfile.aligned.faa
	system("cp $protalnfile $outfile.aligned.faa");
}
else {
	# muscle v3 is the default aligner (program has been tested on muscle v 3.8.31)
	#  -- will use regardless of whether or not --muscle is passed
	system("$muscleexec -in $outfile.unaligned.faa -out $outfile.aligned.faa");
	
	unless ( -e "$outfile.aligned.faa" ) {
		print "Alignment was not generated. Check path to muscle aligner\n\n";
		print "exiting...\n";
		exit;
	}
}

# Convert the aligner output to singleline fasta using awk (see --multiline mode above
# for general awk command). Then mv the singleline fasta to $outfile.aligned.faa
system("awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"%s\",\$0);}  END {printf(\"\\n\");}' < $outfile.aligned.faa | tail -n +2 > temp.faa");
system("mv temp.faa $outfile.aligned.faa");

############################################################################
# Read the amino acid alignment and use to guide the nucleotide alignment
############################################################################
print "\n\nReading aligned amino acid file $outfile.aligned.faa... ";
open (my $PROTF, "$outfile.aligned.faa") or die "Could not open file $outfile.aligned.faa for input.\n";
my @protdata = <$PROTF>; # Read the protein alignment file
close($PROTF) or die "Could not close file $outfile.aligned.faa\n";
my($protlines) = $#protdata + 1;

# Use the number of sites in the protein alignment to calculate the number of codons in
# the codon alignment
chomp($protdata[1]);
my($ncodon) = length($protdata[1]) * 3;

my @protname;
my @protseq;

print "done.\n\nWriting the codon alignments... ";

open (my $AFAF, ">$outfile.codon.alignment.fasta") or die "Could not open file $outfile.codon.alignment.fasta for output.\n";
open (my $PHYF, ">$outfile.codon.alignment.phy") or die "Could not open file $outfile.codon.alignment.phy for output.\n";

# Write the phylip header
print $PHYF "$taxa $ncodon\n";

# Loop over the taxon names in the amino acid alignment
for ($iter=0; $iter<$taxa; $iter++) {

	# Extract the relevant protein sequence information
	$jter = $iter * 2;
	$protname[$iter] = $protdata[$jter];
	chomp($protname[$iter]);
	$protname[$iter] =~ s/>//g;
	$protseq[$iter] = $protdata[$jter+1];
	chomp($protseq[$iter]);
	$len = length($protseq[$iter]);	

	# Loop over the taxon names for the nucleotide data
	for ($kter=0; $kter<$taxa; $kter++) {
		if ( $name[$kter] eq $protname[$iter] ) {
			print $AFAF ">$protname[$iter]\n";
			print $PHYF "$protname[$iter]  ";
			$tval1 = $maxnamelength - length($protname[$iter]);
			for ($mter=0; $mter<$tval1; $mter++) { print $PHYF " "; }
 			$zter = 0;
			for ($lter=0; $lter<$len; $lter++) {
				$aa = substr($protseq[$iter],$lter,1);
				if ( $aa eq "-" ) { 
					print $AFAF "---"; 
					print $PHYF "---"; 
				}
				else {
					$codon = substr($seq[$kter],$zter,3);
					$codon = uc($codon);
					print $AFAF "$codon"; 
					print $PHYF "$codon";
					$zter = $zter + 3;
				} 
			} # end for ($lter=0; $lter<$len; ...
			print $AFAF "\n";
			print $PHYF "\n";
		} # end if ( $name[$kter] eq $protname[$iter] ...
	} # end for ($kter=0; $kter<$taxa ...
} # end for ($iter=0; $iter<$taxa ...

close($AFAF) or die "Could not close file $outfile.codon.alignment.fasta\n";
close($PHYF) or die "Could not close file $outfile.codon.alignment.phy\n";

# Clean up the singleline fasta file if mode is --multiclean or -MC
if ( $multiline == 1 && $clean == 1 ) { system("rm -f $infile"); }

print "done\n\n";
print "Alignments written to:\n";
print "  -- $outfile.codon.alignment.fasta\n";
print "  -- $outfile.codon.alignment.phy\n\n";

print "exiting...\n";


exit;
