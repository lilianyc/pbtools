[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# PBTools
> Protein Block Molecular Dynamics Tools

PBTools is a program to perform analysis on molecular dynamics protein sequences encoded as Protein Blocks. Protein Blocks are structures defined by [De Brevern *et al.*](https://www.ncbi.nlm.nih.gov/pubmed/11025540) to analyze local conformations of proteins.   
The program a simple Python re-implementation of GSATools by [Pandini *et al.*](https://academic.oup.com/bioinformatics/article/29/16/2053/200020) using these Protein Blocks as a structural alphabet.

## Requirements

- Python >= 3.6
- PBXplore
- pandas
- numpy

## Installation

### Developer mode

Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)  

Clone the repository:
```shell
$ git clone https://github.com/lilianyc/projet_court
$ cd projet_court
```

Create conda environment:
```shell
$ conda env create -f environment.yml
```

Remark: for a fully reproducible environment, you could also use:
```shell
$ conda env create -f environment.lock.yml
```
*Conda might give some errors with packaging*

Activate conda environment:
```shell
$ conda activate PBTools
```

Install local package:
```shell
$ pip install -e .
```

## Quickstart

The script can take either [PDB](http://www.wwpdb.org/documentation/file-format) files
or a
[trajectory](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html) and [topology](https://www.mdanalysis.org/docs/documentation_pages/topology/init.html) file.

It will return a csv of an upper triangular Mutual Information matrix for each position of the Protein Block sequence.

### Options

- `-h, --help`:  
Provide a help page

- `-v, --version`:  
Provide the version

- `-o, --output`:  
Give the name of the output file

- `-p, --pdb`:  
Name of PDB file(s) or a directory containing PDB files

- `-x`:  
Name of the trajectory file

- `-g`:  
Name of the topology file

### Examples

After installing PBTools, you can use the Command-Line Interface.

Using
```shell
$ pbtools -o matrix.csv -p 1BTA.PDB
```
will create `matrix.csv`, a file containing the Mutual Information matrix.

## Credits

- GSATools
  - Pandini A, Fornili A, Fraternali F, Kleinjung J. GSATools: analysis of allosteric communication and functional local motions using a structural alphabet. *Bioinformatics* 29(16):2053–2055. https://doi.org/10.1093/bioinformatics/btt326 (2013)

- PBXplore
  - Barnoud J, Santuz H, Craveur P, Joseph AP, Jallu V, de Brevern AG, Poulain P, PBxplore: a tool to analyze local protein structure and deformability with Protein Blocks.
  *PeerJ* 5:e4013 https://doi.org/10.7717/peerj.4013 (2017).

- Protein Blocks
  - A. G. de Brevern, C. Etchebest, S. Hazout,
  Bayesian Probabilistic Approach for Predicting Backbone
Structures in Terms of Protein Blocks. *Proteins* 41:271-87 http://dx.doi.org/10.1002/1097-0134(20001115)41:3%3C271::AID-PROT10%3E3.0.CO;2-Z (2000).
