## 1. 
Documents used for multiple alignment

CLUSTAL 2.1
MUSCLE v3.8.1551 by Robert C. Edgar
MAFFT v7.490 (2021/Oct/30)
Kalign (3.3.1)
T-COFFEE Version_13.41.0.28bdc39 (2019-11-30 10:21:53 - Revision 5d5a1c1 - Build 465)
PRANK v.170427

## 2. 
Code to run 6 possible algorithms
alignments (clustalw,muscle,mafft,kalign,tcoffee,prank) per 10 sequences
DNA (SUP35_10seqs.fa) + parameter variations, if any
were. If everything starts online - links
pages, you can attach screenshots,
if it seems to make sense.

```{bash}
clustalw -INFILE=SUP35_10seqs.fa -OUTPUT=FASTA -OUTFILE=SUP35_10seqs.clustalw.fa
muscle   -in SUP35_10seqs.fa -out SUP35_10seqs_muscle.fa
mafft --auto SUP35_10seqs.fa >SUP35_10seqs_mafft.fa
kalign <SUP35_10seqs.fa >SUP35_10seqs_kalign.fa  #weird format
t_coffee -infile=SUP35_10seqs.fa -outfile=SUP35_10seqs_tcoffee.fa
prank -d=SUP35_10seqs.fa -o=SUP35_10seqs_prank.fa
```

Comparison table over time
work* and comments about the quality
DNA alignment for the above
algorithms**. What is the best algorithm to use?
Necessary
minimum: runtime + alignment length
or a graphical representation (for example,
as in UGENE). Please note that the sequences
code for a protein, it's worth using
when evaluating.

```{bash}
time kalign <SUP35_10seqs.fa >SUP35_10seqs_kalign.fa #time measuring
```

| algorithm  |  time  |  length  |  quality |
| ----------- | ----------- | ----------- | ----------- |
clustalw  |   0m3,052s |    2148bp  |  +++ |
muscle  |   0m2,050s  |   2275bp  |  +++ |
mafft  |  0m1,892s  |   2166bp  |  +++ |
kalign (#weird format)  |  0m0,162s  |  2152  |   + | 
tcoffee (#shows graph in terminal) |  1m28,334s   |  2210bp  |  ++ |
prank   |  Total time 4s  |  2366bp  |  + |

length can be checked in UGENE, as well as the quality. for example, with the ber below the alignment. The grey area in the bar shows % of same nucleatides/aminoacids in each location of alignment. To me they all seem pretty ok, except of prank, and kalign format couldn't be correctly showed in UGENE... So I would choose clustalw/muscle/mafft

## 3. 
SUP35_10seqs_strange_aln.fa is a reverse complementary sequence, so it alignt wrong, though if we run blast on it, everything seems to be fine. To fix that we can use a UGENE function or try these variants i have found:

```{bash}
seqtk seq -r in.fa > out.fa #fastQ
echo ACCTTGAAA | tr ACGTacgt TGCAtgca | rev #without tools
revseq #  EMBOSS
```

## 4. 
Commands / screenshots to run
6 possible alignment options (see
2), but for 250 DNA sequences. Comparative
table with working time and comments
about the quality of the alignment of 250 sequences
DNA (SUP35_250seqs.fa).
Has our choice of algorithm changed?

```{bash}
clustalw -INFILE=SUP35_250seqs.fa OUTPUT=FASTA -OUTFILE=SUP35_250seqs.clustalw.fa
muscle -in SUP35_250seqs.fa -out SUP35_250seqs_muscle.fa
mafft --auto SUP35_250seqs.fa >SUP35_250seqs_mafft.fa
kalign <SUP35_250seqs.fa >SUP35_250seqs_kalign.fa
t_coffee -infile=SUP35_250seqs.fa -outfile=SUP35_250seqs_tcoffee.fa
prank -d=SUP35_250seqs.fa -o=SUP35_250seqs_prank.fa
```

algorithm |  time   | quality |
| ----------- | ----------- | ----------- |
clustalw,   |  ?  |   ?
muscle,  |   53s   |  ++
mafft,  |  0m14,651s  |   +++
kalign,  |  0m2,321s  |  +  |   #weird format
tcoffee,  |  >25m13,742s   |  ? |    #shows graph in terminal
prank  |   72s  |  +

## 5. 
How to add 250 bp to the alignment
two more sequences (SUP35_2addseqs.fsa), previously
aligning them with mafft and muscle?

This trick can be done only in mafft/muscle. At first you should align a new sequense on a one from the made multi alignment, and then insert it in it.

```{bash}
muscle -in SUP35_2addseqs.fa -out SUP35_2addseqs_muscle.fa
muscle -profile -in1 SUP35_250seqs_muscle.fa -in2 SUP35_2addseqs.fa -out SUP35_252seqs_muscle.fa
mafft --auto SUP35_2addseqs.fa > SUP35_2addseqs_mafft.fa
mafft --add SUP35_2addseqs_mafft.fa SUP35_250seqs_mafft.fa > SUP35_252seqs_mafft.fa
```


## 6. 
How to get sequences
amino acids (translate)? Command example
for translation into amino acid sequences.
What problems arise.

```{bash}
transeq -sequence SUP35_10seqs.fa -outseq SUP35_10seqs.t.faa
```

to find an open reading frame (and really working ones): 

```{bash}
getorf -sequence SUP35_10seqs.fa -outseq SUP35_10seqs.g.faa -noreverse -minsize 500
```

bc transeq is stupid, and you might get tottaly wrong proteing if you start from a wrong point in sequence.

## 7. 
Commands / screenshots to run
6 alignment options for
10 protein sequences + variations
parameters, if any.

```{bash}
clustalw -INFILE=SUP35_10seqs.g.faa -OUTFILE=SUP35_10seqs.clustalw.faa -OUTPUT=FASTA -TYPE=protein
clustalo --infile=SUP35_10seqs.g.faa --outfile=SUP35_10seqs.clustalo.faa --verbose
muscle -in SUP35_10seqs.g.faa -out SUP35_10seqs_muscle.faa
mafft --auto SUP35_10seqs.g.faa >SUP35_250seqs_mafft.fa
kalign <SUP35_10seqs.faa >SUP35_10seqs_kalign.faa
t_coffee -infile=SUP35_10seqs.faa -outfile=SUP35_10seqs_tcoffee.faa
prank -d=SUP35_10seqs.faa -o=SUP35_10seqs_prank.faa
```

Comparison table over time
work and comments about the quality
protein alignment. What algorithm
better to use?

algorithm  | time  |  quality
| ----------- | ----------- | ----------- |
clustalw,  |   0m0,329s   |  +++
muscle,  |   0m0,100s  |   ++
mafft,  |  0m0,289s  |   ++
kalign,  |  0m0,161s  |  ?  
tcoffee,  |  1m20,882s   |  ++++  
prank  |   0m3,620s  |  +

## 8. 
Extract from NCBI (using any variations
eutils or screenshots if in a browser)
all on request
18S" (Parapallasea is a taxon and 18S is a gene) and save
in fast file.

```{bash}
esearch -db nucleotide -query "Parapallasea 18S" | efetch -format fasta >Parapallasea_18.fa
muscle -in Parapallasea_18.fa -out Parapallasea_18.fa.muscle.aln
mafft --auto Parapallasea_18.fa > Parapallasea_18.fa.mafft.aln
```

## 9.
What goes wrong with sequence alignment
in the file Parapallasea_18S.fa and with what parameters you can
get the right answer?

The obtained sequences were partial, and they were from different parts of the gene. Thus, we have to merge them into one sequence with cap3, and then it should be aligned. Also it can be done with EMBOSS or UGENE (consensus -> export). Without that the alignment looks weird, but it is possible to understand it more or less...

## 10. 
Commands for shaping
from the sequence set Ommatogammarus_flavus_transcriptome_assembly.fa
base for blast, and for searching in this base
protein sequence of Acanthogammarus_victorii_COI.faa with
recording the results in a table (text
tab delimited). Attention: origin
mitochondrial sequences.
What is important to consider when searching?

Extract sequence with best
match in a separate file.

```{bash}
makeblastdb -in Ommatogammarus_flavus_transcriptome_assembly.fa -dbtype nucl -parse_seqids
tblastn -query Acanthogammarus_victorii_COI.faa -db Ommatogammarus_flavus_transcriptome_assembly.fa -outfmt 6 #for mitochondrial genes
## fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
blastdbcmd -db Ommatogammarus_flavus_transcriptome_assembly.fa -entry TRINITY_DN8878_c0_g1_i2 -out Ommatogammarus_flavus_COI.fa #writing a sequence to file
```
