codonpair
=========

`codonpair` calculates codon pair score and codon pair bias. CPS
values are close, but not identical, with those produced by the perl
script from Dimitris Papamichail (`cps_perl` directory) and, presumably,
used in the following work: 

Virus attenuation by genome-scale changes in codon pair bias.
Coleman JR1, Papamichail D, Skiena S, Futcher B, Wimmer E, Mueller S.
Science. 2008 Jun 27;320(5884):1784-7. doi: 10.1126/science.1155761.
https://www.ncbi.nlm.nih.gov/pubmed/18583614

The difference(s) between the (python) implementation here and the earlier
work requires further investigation. 2 regression tests don't pass.

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

