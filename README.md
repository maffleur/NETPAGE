# NETPAGE
NETPAGE (NETwork Propagation-based Assessment of Genetic Events) is a novel computational framework for gene-based association testing of rare variants that integrates prior knowledge about tissue-specific gene interaction networks (Scelsi et al., PLOS Comp Biol 2021, in press; DOI: 10.1371/journal.pcbi.1008517).

The framework combines network propagation with sparse regularised regression. Here we provide a Python module to perform network propagation. 

## Getting started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
You will need Python 2.7 to use this module. A list of module dependencies can be found in the `requirements.txt` file.

### Installation
No installation is required. Just clone this repo to your machine, make sure you have all the packages needed, and you're ready to go!

### An example
Try reproduce the results in the `example` folder by running:
```
python RunNetworkSmoothingParser.py -i example/example-net -g example/example-geneburden.txt --typ burden --alpha 0.5
```
and
```
python RunNetworkSmoothingParser.py -i example/example-net -g example/example-geneburden.txt --typ burden --alpha 0.9
```
*NOTE:* Please make sure that your gene burden file contains a column named "PTID" containing the subject IDs.
### Help page
```
python RunNetworkSmoothingParser.py -h

```
