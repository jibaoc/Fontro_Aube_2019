# Description

Fasta generator program, allow to create a fasta file (containing DNA sequences) with your desired number of sequences.
You can specify the proportion of alanine, guanine, thymine and cytosine in your sequences

# Prerequisite

This program works with `python 2.7` and the following modules :
* `random` : for random suffling of sequences
* `argparse` : for parsing arguments
* `os` : to check if a directory exists
* `math` : to truncate numbers

# Running the program

A description of all possible parameters can be displayed by running this command :

```sh
python2 src/fasta_generator.py --help
```

## Example of command

```sh
# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of alanine equals to 80%
# and a proportion to G, T, C equal to 6.67 %
python2 src/fasta_generator.py --prop_A 0.8 --filename "my_fasta" --nbr_seq
100 --size_inf 50 --size_max 300

# Creation of a fasta file my_fasta.fasta in your current working directory
# It will contain 100 sequences having a size between 50 and 300 nt.
# Every sequences will have a proportion of alanine equals to 50%
# and a proportion to G and  T equal to 20 %
# and a proportion of C of 10%
python2 src/fasta_generator.py --prop_A 0.5 --prop_C 0.1 --prop_T 0.2
--filename "my_fasta" --nbr_seq 100 --size_inf 50 --size_max 300

```
