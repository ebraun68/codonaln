# codonaln
Programs to align sequences while respecting codon structure.

Two perl programs (one to check frames and one to perform the alignment) and two
awk programs to convert file formats.

--------------------------------------------------------------------------------
# checkframe.pl

A simple perl program to examine a set of sequences and identify the best
reading frame. Reads a table of nucleotide sequences and count the number of
stop codons. Outputs a fasta file with each sequence "padded" with N's to put the
sequence in frame (with frame information appended to the sequence name) and a
tab-delimited file with the numbers of stop codons in each frame as well as the
number of "problematic" codons.

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

