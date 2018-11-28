# ISAErrorCorrect
Read, characterize, combine ISA alignments

This program is the fourth step in the ISA pipeline.

1: ISAFil1 filters fastsq/fastq reads and combines identical sequences -> .fasta
2: BLAT alignes reads to target genome -> .blast8
3: blast8toAln.py filters BLAT output -> .aln
*4: ISAErrorCorrect reads in .aln and .fasta. Combines reads with matching alignments and does error correction -> .tran
5: buildGlobalISA reads .tran files for all tests of the same subject. 
  Does some further error correction and normalizes read counts of all samples
