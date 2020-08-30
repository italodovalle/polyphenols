# Network Medicine Framework to Predict Health Impact of Dietary Polyphenols

Jupyter (IPython) Notebook and required files for the proximity-based analysis in the *"Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework"* manuscript.

## Instructions

```
git clone https://github.com/italodovalle/polyphenols.git
cd polyphenols
git init submodule
git submodule update
```


## Usage

**run_proximity.py**

```
usage: run_proximity.py [-h] -i INTERACTOME -d DISEASEGENES -t TARGETS
                        [-c CORES] [-r NRANDOM] -o OUTFOLDER [-sn SOURCE_NODE]
                        [-tn TARGET_NODE] [-lcc LCC] [-sep SEPARATOR]
                        [-test TEST_RUN]

Run Proximity

optional arguments:
  -h, --help       show this help message and exit
  -i INTERACTOME   interactome file
  -d DISEASEGENES  disease genes file
  -t TARGETS       target genes file
  -c CORES         number of cores to use, default = 6
  -r NRANDOM       n random permutations, default = 1000
  -o OUTFOLDER     outfolder
  -sn SOURCE_NODE  name source node, default = proteinA
  -tn TARGET_NODE  name target node, default = proteinB
  -lcc LCC         LCC, default = True
  -sep SEPARATOR   separator, default = ","
  -test TEST_RUN   test run, default = False
```


## Data Files

The data folder contains:

* data/HumanInteractome_v2017.txt: Human interactome (see manuscript for details)
* data/GenesDisease.tsv: Disease-gene associations for MeSH disease terms curated in Menche et al. (2015, Science)
* data/PolyphenolProteinInteractions.csv: Polyphenol protein targets obtained from the STITCH database
* data/connectivity_map: Perturbation profiles obtained in the Connectivity Map database (https://clue.io/) for selected polyphenols

## Citation

do Valle, I. F. et al. Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework. bioRxiv 2020.08.27.270173 (2020) doi:10.1101/2020.08.27.270173.
