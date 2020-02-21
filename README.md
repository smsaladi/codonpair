[![Build Status](https://travis-ci.org/smsaladi/codonpair.svg?branch=master)](https://travis-ci.org/smsaladi/codonpair)
[![PyPI version](https://badge.fury.io/py/codonpair.svg)](https://badge.fury.io/py/codonpair)
![PyPI - Downloads](https://img.shields.io/pypi/dm/codonpair)
[![DOI](https://data.caltech.edu/badge/77872126.svg)](https://data.caltech.edu/badge/latestdoi/77872126)

codonpair
=========

`codonpair` calculates codon pair score and codon pair bias. CPS
values are identical to those produced by the perl script from
Dimitris Papamichail (`cps_perl` directory) and, presumably,
used in the following work: 

    Virus attenuation by genome-scale changes in codon pair bias.
    Coleman JR1, Papamichail D, Skiena S, Futcher B, Wimmer E, Mueller S.
    Science. 2008 Jun 27;320(5884):1784-7. doi: 10.1126/science.1155761.
    https://www.ncbi.nlm.nih.gov/pubmed/18583614


## Installation

Either, clone the repo and install with pip

```shell
git clone git@github.com:smsaladi/codonpair.git
pip install ./codonpair
```

Or... have pip handle the details:

```shell
pip install git+git://github.com/smsaladi/codonpair@master#codonpair
```

All dependencies should be checked for and, if necessary, installed
automatically by `pip`.


## Usage

Initialize a `codonpair.CodonPair` object by specifying a list of reference sequences
`CodonPair.from_sequences`, from a named reference `CodonPair.from_named_reference`,
a reference file `CodonPair.from_reference_file`,
or simply providing a `pd.DataFrame` with codon counts to `CodonPair`.

The following named references are recognized/bundled with this package.

* `E. coli` (BL21 DE3)
* `S. pneumoniae` (TIGR4)
* `cps_perl` -  the reference file provided with the perl implementation

The default constructor `CodonPair()` uses the `E. coli`.


Then calculate the codon pair score for a provided sequence with `CodonPair.cpb`
which returns a dictionary with the

* total codon pair score `total_cps` - the sum of the values of each codon pair
* the number of codons `n_pair` - excluding codon pairs not found in the reference
* the codon pair bias `cpb` - `total_cps/n_pair`

For one-off calculations, `codonpair.calc_cpb` can be used directly
for with the sequence of interest (calling the default constructor under the hood).

```python
import codonpair
cp = codonpair.CodonPair.from_named_reference('E. coli')
cp.cpb("ATGATCCCCTTACAACATGGACTGATCCTCGCGGCAATCTTATTCGTTCTTGGCTTAACC")
```

For convenience, the executable `cps` installed into the path by pip:

```bash
cps test.fasta > test.scores.txt
```

See `CodonPair.write_reference` to write codon pair counts for a reference set to
the filename provided to be used with future calculations.
