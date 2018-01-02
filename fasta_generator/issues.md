<div class="alert alert-info">

**:warning:** There are serveral ***issues*** about the random generation procedure described in *command.md* file :
1. The generation of random fasta sequences with di-nt frequency only control the frequencies of half dnt in the sequences :
  If *seq=AC and next dnt is GT the sequence will become ACGT and the dnt CG wans't set by the program !*
2. The creation of random fasta with a codon proportion that match CCE/ACE/ALL exons doest not control the possible bias at hexanucleotides level
3.  The creation of random fasta with a codon proportion that match CCE/ACE/ALL with a given di-nt frequency can oly create sequence having this dn-nt frequnecy greater than coincidence. You can't use this to impoverished a sequence of some di-nucleotides
</div>
