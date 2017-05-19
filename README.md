# disseqd
Detection of Interspersed Signatures in SEQuence Data - A software to detect subsequences within sequences, based on distinct kmer frequencies.

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [disseqd](#disseqd)
	- [Introduction](#introduction)
	- [Documentation](#documentation)
	- [Setup](#setup)
	- [Usage](#usage)
	- [Tips and tricks](#tips-and-tricks)
	- [Questions and issues](#questions-and-issues)
	- [Contributing](#contributing)
	- [Release history](#release-history)
	- [Team](#team)
	- [Credit](#credit)
	- [License](#license)
	- [Links](#links)

<!-- /TOC -->

## Introduction
`disseqd` implements a hidden Markov model to label kmers within a sequence with distinct sequence categories. The idea is to identify subsequences of different origin or with specific properties, which can be distinguished by their kmer composition.

## Documentation
Input and output files of `disseqd` are formatted in fasta-style:
```python
>option
value
[value]
[...]
>option
value
[value]
[...]
```
### Input
#### Data
`disseqd` can be used to disseqd the observed data, or to train the model parameters from additional input sequences prior to disseqdion. For both purposes, `disseqd` requires a data file to be disseqded:
```python
>observed_data
GCCGTCACCGCCGGCCGTCCCGCGCCCTGCCGAGAGCCACACCCTGGCCTGTCACCATACCCATCC
```
#### Model
`disseqd` also needs the model parameters to be used. These are simplest specified via a model file, but (their short versions) can also be supplied via the command line with the corresponding flags (see `python disseqd -h`), and do then override that option if already set in the model file.

Two single-choice options can be specified:
```python
>kmer # number of consecutive characters to consider
2
>alphabet # legal characters
ACGT
```

And four multi-choice options can be specified. Short version:
```python
>nstates # number of states
2
>start # start probability for all states
0.5
>transitions # expected lengths of kmer tracts per state
10
20
>emissions # sequence files for obtaining emission probabilities
examples/cpg.txt
examples/noncpg.txt
```
These short options can also be specified on the command line with the corresponding flag (see `python disseqd -h`), and do then override that option if already set in the model file.

Long version:
```python
>nstates # names of states
CPG
NCPG
>start # start probabilities for each state
0.4
0.6
>transitions # transition probabilities between states
0.9	0.1
0.05	0.95
>emissions # emission probabilites of kmers per state
AA	0.00025	0.00025
AC	0.99925	0.49975
AG	0.00025	0.25
AT	0.00025	0.25
CA	0.0588235294118	0.33325
CC	0.529411764706	0.53305
CG	0.352941176471	0.00025
CT	0.0588235294118	0.13345
GA	0.00025	0.333333333333
GC	0.66625	0.333333333333
GG	0.11125	0.166666666667
GT	0.22225	0.166666666667
TA	0.00025	0.20005
TC	0.66625	0.39985
TG	0.33325	0.39985
TT	0.00025	0.00025
```
The model file can contain any of the 6 different parameters specified above. Aside from kmer and alphabet, the short version of the model file can contain the number of states, the start probability for all states, the expected lengths of kmer tracts per state, and sequence files for obtaining emission probabilities. The latter four options can also be specified in long format, and can then contain the names of states, the start probabilities for each state, the transition probabilities between states, and the emission probabilities of kmers per state.

### Output
#### Decoding
PREFIX.decode.txt contains the decoded observations:
```python
>loglikelihood:-72.21101909	CPG:0	NCPG:1
000000000000000000000000000000011111111111111111111111111111111111
```
#### Model
PREFIX.model.txt contains the full model, long version (see Input).

## Setup
Running the software requires python 3 installed on your system.

To get `disseqd` from github, type the following at a command line:
```
git clone -b develop https://github.com/charlgren/disseqd.git
```

## Usage
After acquiring the software, `disseqd` needs to be called from a terminal, with input specified as in the Documentation section.

Basic usage:
```
python disseqd.py -d data.fa -m model.txt
```
For more options:
```
python disseqd.py -h
```

## Tips and tricks


## Questions and issues


## Contributing


## Release history


## Team


## Credit


## License
Released under the [MIT License](https://github.com/charlgren/disseqd/blob/develop/LICENSE).

## Links
