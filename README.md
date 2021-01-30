# Network Medicine Framework to Predict Health Impact of Dietary Polyphenols

Jupyter (IPython) Notebook and required files for the proximity-based analysis in the *"Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework"* manuscript.


## Table of Contents

  * [Instructions](#instructions)
  * [Network Proximity](#network-proximity)
  * [Implicit Associations](#implicit-associations)
  * [Gene Set Enrichment Analysis](#gene-set-enrichment-analysis)
  * [Data Files](#data-files)
  * [Citation](#citation)


## Instructions

**Clone repository**

```
git clone https://github.com/italodovalle/polyphenols.git
```


**Install requirements**

```
cd polyphenols/
pip install -r requirements.txt
```


## Network Proximity

**run_proximity.py**

```
usage: run_proximity.py [-h] -i INTERACTOME -d DISEASEGENES -t TARGETS
                        [-c CORES] [-r NRANDOM] [-o OUTFOLDER]
                        [-sn SOURCE_NODE] [-tn TARGET_NODE] [-lcc LCC]
                        [-sep SEPARATOR] [-test TEST_RUN] [-chem CHEMICAL_COL]
                        [-prot PROTEIN_COL] [-tmp TMP_DIR]

Calculate Network Proximity

optional arguments:
  -h, --help          show this help message and exit
  -i INTERACTOME      interactome file
  -d DISEASEGENES     disease genes file
  -t TARGETS          chemical protein interactions
  -c CORES            number of cores to use, default = 6
  -r NRANDOM          n random permutations, default = 1000
  -o OUTFOLDER        outfolder
  -sn SOURCE_NODE     name source node Interactome, default = proteinA
  -tn TARGET_NODE     name target node Interactome, default = proteinB
  -lcc LCC            Consider only Interactome LCC, default = True
  -sep SEPARATOR      separator Interactome, default = ","
  -test TEST_RUN      test run
  -chem CHEMICAL_COL  chemical column in chemical-protein interaction file
  -prot PROTEIN_COL   chemical column in chemical-protein interaction file
  -tmp TMP_DIR        temporary folder
```

### Examples

**Example**

`./run_proximity.py -i data/HumanInteractome_v2017.csv -d data/PeroxisomalDisorders.csv -t data/Genistein.csv -r 1000 -sn EntrezA -tn EntrezB -chem chemical -prot entrez_id`

* run time: 30 seconds


**Entire dataset**

`./run_proximity.py -i data/HumanInteractome_v2017.csv -d data/GenesDisease.csv -t data/PolyphenolProteinInteractions.csv -c 6 -r 1000 -tmp tmp/ -sn EntrezA -tn EntrezB -test True -chem chemical -prot entrez_id`

* run time: it can take several days to calculate the proximity for the entire dataset.


## Implicit Associations

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
python expand_associations.py -m data/mtrees2018.bin -i data/ctd_disease_chemical_phenolexplorer_therapeutic.csv -o data/ctd_polyphenols_implicit_explicit.csv
```

## Gene Set Enrichment Analysis


**run_gsea.py**

Usage


```
usage: run_gsea.py [-h] -indir INDIR -i INTERACTOME -d DISEASEGENES -t TARGETS
                   [-c CORES] [-r NRANDOM] [-o OUTFOLDER] [-sn SOURCE_NODE]
                   [-tn TARGET_NODE] [-lcc LCC] [-sep SEPARATOR]
                   [-chem CHEMICAL_COL] [-prot PROTEIN_COL] [-tmp TMP_DIR]
                   [-cell CELL_LINE]

Calculate Enrichment Score

optional arguments:
  -h, --help          show this help message and exit
  -indir INDIR        folder with .gct files
  -i INTERACTOME      interactome file
  -d DISEASEGENES     disease genes file
  -t TARGETS          chemical protein interactions
  -c CORES            number of cores to use, default = 6
  -r NRANDOM          n random permutations, default = 1000
  -o OUTFOLDER        outfolder
  -sn SOURCE_NODE     name source node Interactome, default = proteinA
  -tn TARGET_NODE     name target node Interactome, default = proteinB
  -lcc LCC            Consider only Interactome LCC, default = True
  -sep SEPARATOR      separator Interactome, default = ","
  -chem CHEMICAL_COL  chemical column in chemical-protein interaction file
  -prot PROTEIN_COL   chemical column in chemical-protein interaction file
  -tmp TMP_DIR        temporary folder
  -cell CELL_LINE     restrict to cell line, default = MCF7

```


Example

```
./run_gsea.py -indir data/connectivity_map/ -i data/HumanInteractome_v2017.csv -d data/PeroxisomalDisorders.csv -t data/Genistein.csv -r 10 -sn EntrezA -tn EntrezB -chem chemical -prot entrez_id
```

* Results written in `out_gsea`
* Run time: It take hours/days for the entire dataset

## Data Files

The data folder contains:

* `data/HumanInteractome_v2017.csv`: Human interactome (see manuscript for details)
* `data/GenesDisease.csv`: Disease-gene associations for MeSH disease terms curated in Menche et al. (2015, Science)
* `data/PolyphenolProteinInteractions.csv`: Polyphenol protein targets obtained from the STITCH database
* `data/connectivity_map`: Perturbation profiles obtained in the Connectivity Map database (https://clue.io/) for selected polyphenols
* `data/mtrees2018.bin`: MeSH hierarchy
* `ctd_disease_chemical_phenolexplorer_therapeutic.csv`: all disease-polyphenol associations obtained by combining data from the databases CTD and PhenolExplorer
* `ctd_polyphenols_implicit_explicit.csv`: output file of the code `expand_associations.py`
* `Genistein.csv`: Genistein targets (for tests)
* `PeroxisomalDisorders.csv`: Peroxisomal Disorders disease genes (for tests)

## Citation

do Valle, I. F. et al. Predicting the Health Impact of Dietary Polyphenols Using a Network Medicine Framework. bioRxiv 2020.08.27.270173 (2020) doi:10.1101/2020.08.27.270173.
