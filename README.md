# codonaln
Programs to align sequences while respecting codon structure.

Two perl programs (to perform the alignment and one to check frames) and two
awk programs to convert file formats. Both should work in *nix environments, 
including MacOS.

--------------------------------------------------------------------------------
### codon_align.pl

A perl program that translates a set of nucleotide sequences passed as a fasta file,
aligns the amino acid sequences using muscle (https://www.drive5.com/muscle/) or mafft 
(https://mafft.cbrc.jp/), and uses the amino acid alignment to produce a codon level
alignment. The program defaults to version 3 of muscle (although it can also be used
with version 5 of muscle) The path to the executables of the alignment programs must be 
muscle or mafft, although this can be changed (see below). The program can also use a
protein sequence alignment supplied by the user.

The philosopy of this program is to tolerate ambiguous codons and stop codons. They 
are translated as X in the file used as input for the aligner. This philosophy means
that the user must check the reading frames of the input file of unaligned nucleotide
sequences. The program assumes the first position of the input nucleotide fasta file
is a first codon position. codon_align.pl does ot assume that the lengths of the input 
sequences are a multiple of three; it will automatically pad the 3' end of any sequences
that are not a multiple of three (this is likely to create an ambiguous codon, but it 
allows any sequences to be aligned).

checkframe.pl can can be used to check sequences and it will output a file with the ends 
of the sequences padded with N's to create seqeunces with the correct lengths (checkframe.pl 
chooses the appropriate frame by minimizing the number of in frame stop codons). Note that
the tolerance of stop codons means that the alignments are useful for phylogenetic analyses 
of nucleotides but some analytical programs (e.g., programs to estimate the dN/dS ratio like
codeml in PAML) do not tolerate in frame stop codons. Therefore, users should exercise 
caution if there are any stop codons present in the alignment. 

codon_align.pl uses the standard genetic code: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1), 
although an alternative genetic code can be substituted by altering the hash called %stdgencode. 
checkframe.pl focuses on stop codons and the relevant set of stop codons is supplied by 
the user, so it does not use a specific code.

Minimal usage of codon_align.pl is:

```
perl codon_align.pl input.fasta output
```

The minimal command defaults to muscle v3 as an aligner and it assumes that the fasta input 
file has a single line, like this:

```
>Example_sequence
NTGNTTTNTTTTCATTTCAGCAGAGTATTTTTTNTCACCATCATGGCCTATGACCGGTACATTGCCATCTGCAAACC
```

If your input fasta file has endlines embedded like this:

```
>Example_sequence
NTGNTTTNTTTTCATTTCAGCAGAGTATTTTTTNTCACCATCATGGCCTA
TGACCGGTACATTGCCATCTGCAAACC
```

You should use the multiline fasta mode or clean multiline fasta mode:

```
perl codon_align.pl input.fasta output -M
```
or
```
perl codon_align.pl input.fasta output -MC
```

Both multiline fasta modes will convert the input fasta file to a single line fasta file
before proceeding with the translation and alignment. The clean multiline fasta mode will
remove the single line fasta file. Note that multiline mode does not alter a single line 
fasta file, so routinely using -MC will allow the program to use either type of file without 
saving any new files.

#### Output files

The program has five output files. They are named using the <outfile> name passed to the 
program and different extensions:

```
<outfile>.unaligned.faa         -- amino acid sequences (fasta format)
<outfile>.aligned.faa           -- aligned amino acid sequences (fasta format)
<outfile>.seq_summary.txt       -- summary of sequence lengths and errors
<outfile>.codon.alignment.fasta -- codon aligned sequences (fasta format)
<outfile>.codon.alignment.phy   -- codon aligned sequences (phylip format)
```

Note that the phylip format codon alignment is actually relaxed phylip format (long sequence names 
are tolerated) with all nucleotides on a single line. 

The seq_summary.txt file is a tab-delimited text file with five columns: 1. Sequence name (Name), 
2. the number of nucleotides in the sequence (sites), 3. whether the input sequence was a multiple 
of three nucleotides (EndFrame), 4. number of in frame stop codons (Nstops), and 5. number of 
ambiguous codons and stop codons (NumX).


#### Alignment programs

The default alignment program is version 3 of muscle (the program has been tested with 
muscle v3.8.31). The system call for muscle is: muscle -in unaligned.faa -out aligned.faa 
(where unaligned.faa are aligned.faa the names for the protein fasta files). The other 
alignment programs are used as follows. 

Note that the paths to the aligners (stored in the variables $muscleexec, $muscle5exec, and 
$mafftexec) may need to be changed, depending on the locations and names of the relevant
programs on your system. Also note that I have assumed that the program is run in -MC mode, 
which is robust to single line or multiline fasta files.

#### 1. mafft -
Align using the command: mafft --auto unaligned.faa > aligned.faa

```
perl codon_align.pl input.fasta output --mafft -MC
```

#### 2. muscle v5 with align mode -
Align using the command: muscle -align unaligned.faa -output aligned.faa

```
perl codon_align.pl input.fasta output --muscle5 -MC
```

#### 3. muscle v5 with super5 mode -
Align using the command: muscle -super5 unaligned.faa -output aligned.faa
(-super5 can be used with large alignments when -align is too slow)

```
perl codon_align.pl input.fasta output --super5 -MC
```

#### 4. User supplied protein alignment -

The aligner can also be used with user supplied protein sequence alignment, as follows:

```
perl codon_align.pl input.fasta output --protaln TESTatp1.protaln.faa -MC
```

The protein alignment passed with --protaln must be in aligned fasta format and it must 
not have the same name as the the aligned protein outfile (the program will print an error 
message and exit if this is the case). The program will check whether the alignment file 
exists but it will not perform more detailed checks. If the codon aligment appears problematic 
you should check the protein alignment passed with --protaln.

#### Help message

If you call the program without the minimal command line arguments the program will print 
a help message and exit:

```
Usage:
  $ codon_align.pl <infile> <outfile> OPTIONAL: [aligner] [--protaln filename] [-M or --multiline]
   The first two arguments are mandatory
      infile  = single line aligned fasta infile
      outfile = prefix for output files
      - This program will write five output files:
          <outfile>.unaligned.faa         -- amino acid sequences
          <outfile>.aligned.faa           -- aligned amino acid sequences
          <outfile>.seq_summary.txt       -- summary of sequence lengths and errors
          <outfile>.codon.alignment.fasta -- codon aligned sequences (fasta format)
          <outfile>.codon.alignment.phy   -- codon aligned sequences (phylip format)
   The third and fourth arguments are optional
      - muscle (v3) is the default aligner (can also specified using --muscle3)
      - mafft can be specified using --mafft
      - muscle v5 align mode can be specified using --muscle5
      - muscle v5 super5 mode can be specified using --super5
      - pass --protaln followed by a filename to use an existing protein alignment 

      - use --translate to translate sequences without aligning
      - use --check to assess numbers of in frame stop codons without aligning sequences
      - use --translate to translate sequences without aligning

      - If -M or --multiline is passed the input file can be a multiline fasta file
      - The multiline fasta file will be converted to singleline and saved with the
        the filename infile.singleline
      - If you want to use a multiline input file but do not want to save a the singleline
        file use the clean multiline mode by passing -MC or --multiclean
exiting...
```

--------------------------------------------------------------------------------
### checkframe.pl

A simple perl program to examine a set of sequences and identify the best reading
frame for each sequence. The best frame is chosen by minimizing the number of stop
codons.

This program reads a tab-delimited table of nucleotide sequences and counts the number
of stop codons. The tabular input format can be generated from a fasta file using the 
fasta2tbl.awk program (see below). The program outputs a fasta file with the best frame
for each sequence indicated. If necessary, the ends of the sequences in the fasta output
file will be "padded" with N's to create a sequence where the first nucleotide is the
first position of a codon and the sequence lenght is a multiple of three. The program
also generate a tab-delimited file with the numbers of stop codons in each frame as well 
as the number of "problematic" codons (see below for details).

Typical usage is:

```
perl checkframe.pl input.tbl output STD
```

This reads the input file (input.tbl), counts the number of stop codons in each frame
(assuming the standard genetic code), and outputs a fasta file (output.bestframe.fasta)
with each sequence in the best frame as well as a tab-delimited file (output.frames.txt)
with the numbers of stop codons in each frame. The "best frame" is the one with the 
minimum number of stop codons. Full usage is:

```
  $ checkframe.pl <infile> <outfile> <STOP1> <STOP2> <STOP3> <STOP4>
   The first two arguments are mandatory
      infile  = tab-delimited table of sequence names and sequences
      outfile = prefix for output files
      - This program will write two output files:
          <outfile>.bestframe.fasta  -- nucleotide sequences in correct frame
          <outfile>.frames.txt       -- numbers of stop codons in each frame
      STOPx   = list of stop codons to consider
      - Use STANDARD or STD for the standard genetic code (TAA, TAG, TGA)
      - Use VMT for the vertebrate mitochondrial code (TAA, TAG, AGA, AGG)
      - Otherwise list the stop codons to consider
          For example:
            TAA, TAG, TGA -- Standard code without ambiguities
            TGA           -- Ciliate, Dasycladacean, and Hexamita nuclear code
          Up to eight stop codons (including ambiguities) can be specified
```

The input format is a table with two columns. The first column is the taxon name
(without any internal whitespace) and the second is the nucleotide sequence. A fasta
file can be converted into this tabular format as follows:

```
awk -f fasta2tbl.awk fastafile > tblfile
```

The .bestframe.fasta output file has the best frame added to each name and, if the best
frame is not in the first forward frame, the 5' end of the sequence will be padded with N's.
If the sequence length is not a multiple of 3 the 3' end of the sequence will be padded
with N's.

The .frames.txt output file is tab-delimited text. Column 1 is the sequence name, columns
2-7 are the numbers of problem codons for the 3 forward frames, columns 8-13 are the 
numbers of problem codons for the 3 reverse frames, Column 14 is the best frame (the frame 
with the minimum number of stop codons) and column 15 is the minimum number of stop codons 
in any frame. "stop" is the number of stop codons (ideally 0 or 1) and "gap" is the number 
codons with one or two gaps (the "gap" counts is useful for assessing alignments it will be 
0 for unaligned sequences).

In full, the 15 columns are:
1.  Name
2.  F1stop
3.  F1gap
4.  F2stop
5.  F2gap
6.  F3stop
7.  F3gap
8.  R1stop
9.  R1gap
10. R2stop
11. R2gap
12. R3stop
13. R3gap
14. Best_frame
15. MinStops

The input file can also be a single line relaxed phylip file without the header (i.e.,
the file produced by tail -n +2 input.phy > input.tbl).

