# Network Medicine Framework to Predict Health Impact of Dietary Polyphenols

Jupyter (IPython) Notebook and required files for the proximity-based analysis in the *"Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework"* manuscript.

## Instructions

```
git clone https://github.com/italodovalle/polyphenols.git
cd polyphenols
git init submodule
git submodule update
pip install -r requirements.txt
```


## Network Proximity

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

## Obtaining implicit disease-chemical associations associations

* Consider the example of the disease 'Cardiovascular diseases' (C14.XXX):
    * The following chemicals have the 'therapeutic' label in CTD: Genistein, Kaempferol, Resveratrol, Daidzein
* Now, if you look at the specific diseases that are under 'Cardiovascular diseases' in the MeSH hierarchy we find the following chemicals with 'therapeutic' labels in CTD:
    * cardiomyopathies C14.280.238 (Narigenin, Apigenin)
    * cardiomyopathy, hypertrophic C14.280.484.150.070.160 (Quercetin)
    * heart diseases C14.280 (-)-Epigallocatechin 3-O-gallate
    * etc
* Therefore, we can consider, for example, that the association 'narigenin' and 'cardiovascular diseases' is implicit
* Therefore, we mapped all explicit associations in the MeSH hierarchy (4816 diseases) to expand explicit to implicit associations
* MeSH mapping to the different branches:
    * all branches were considered
* Considering only diseases under the branch C (Diseases). Some phenotypes in CTD might be under the branch F (Psycology and Psychiatry)
* MeSH Supplementary Concept Data ID were manually mapped
* Input
    * `data/mtrees2018.bin`
    * `data/ctd_disease_chemical_phenolexplorer_therapeutic.csv`
* Output:
    * `data/ctd_polyphenols_implicit_explicit.csv`


**expand_associations.py**

Usage

```
usage: expand_associations.py [-h] -m MESH -i INFILE [-o OUTFILE]

Expand associations in MeSH tree

optional arguments:
  -h, --help  show this help message and exit
  -m MESH     MeSH File
  -i INFILE   association table
  -o OUTFILE  outfile
```


Example:

```
python expand_associations.py -m ../data/mtrees2018.bin -i ../data/ctd_disease_chemical_phenolexplorer_therapeutic.csv -o ../data/ctd_polyphenols_implicit_explicit.csv
```


## Data Files

The data folder contains:

* data/HumanInteractome_v2017.txt: Human interactome (see manuscript for details)
* data/GenesDisease.tsv: Disease-gene associations for MeSH disease terms curated in Menche et al. (2015, Science)
* data/PolyphenolProteinInteractions.csv: Polyphenol protein targets obtained from the STITCH database
* data/connectivity_map: Perturbation profiles obtained in the Connectivity Map database (https://clue.io/) for selected polyphenols
* data/mtrees2018.bin: MeSH hierarchy
* ctd_disease_chemical_phenolexplorer_therapeutic.csv: all disease-polyphenol associations obtained by combining data from the databases CTD and PhenolExplorer

## Citation

do Valle, I. F. et al. Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework. bioRxiv 2020.08.27.270173 (2020) doi:10.1101/2020.08.27.270173.
