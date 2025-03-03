> Laboratory TIMC / MAGe Grenoble

# GBA centrality

This repository contains scripts for _GBA (Guilt-by-association) centrality_.

_GBA (Guilt-by-association) centrality_ is a network propagations algorithm for disease gene prioritization. The method assigns scores to genes that represent their likelihood of being causal for the phenotype. It takes into account the topology of the protein-protein interaction network (interactome) and prior knowledge about causal genes for the phenotype of interest.

## How to construct the human interactome

### Step 1. Download interactome data

Download Uniprot data

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

Download PPI data from BioGRID

```
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
```

```
unzip BIOGRID-ORGANISM-LATEST.mitab
```

Download PPI data from IntAct

```
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
```

```
unzip intact.zip
```

Download PPI data from Reactome

```
wget https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.psi-mitab.txt
```


### Step 2. Parse interactome data

Parse Uniprot

```
gunzip -c uniprot_sprot.dat.gz | python Interactome/uniprot_parser.py > uniprot_parsed.tsv
```

Parse BioGRID

```
python Interactome/interaction_parser.py --interaction_file BIOGRID-ORGANISM-Homo_sapiens*.mitab.txt --uniprot_file uniprot_parsed.tsv > interactions_Biogrid.tsv
```

Parse IntAct

```
python Interactome/interaction_parser.py --interaction_file intact.txt --uniprot_file uniprot_parsed.tsv > interactions_Intact.tsv
```

Parse Reactome

```
python Interactome/interaction_parser.py --interaction_file reactome.homo_sapiens.interactions.psi-mitab.txt --uniprot_file uniprot_parsed.tsv > interactions_Reactome.tsv
```


### Step 3. Build the high-confidence human interactome

```
python Interactome/build_interactome.py \
  --interactions_parsed_files interactions_Biogrid.tsv interactions_Intact.tsv interactions_Reactome.tsv \
  --uniprot_file uniprot_parsed.tsv > interactome_human.sif
```


### Step 4. Prepare a list of causal genes

Create a TSV file `causalGenes.tsv` (without a header) with 2 columns: gene name, pathology


## How to run _GBA centrality_

As input, GBA centrality takes an interactome SIF file, a Uniprot file and a TSV file with known causal genes for the phenotype of interest.

Example usage for an infetility phenotype MMAF and parameters alpha=0.5, d_max=10

```
python GBA_centrality.py \
  -i interactome_human.sif \
  --uniprot_file uniprot_parsed.tsv \
  --causal_genes_file causal_genes_infertility.tsv \
  --alpha 0.5 \
  --d_max 10 \
  --patho MMAF \
  1> output/scores.tsv \
  2> output/log.txt
```

### Python environment

GBA centrality is written in Python :snake: and requires the following dependencies: [NumPy](https://numpy.org/) and [Networkx](https://networkx.org/).

We recommend installing them via [Python venv](https://docs.python.org/3/library/venv.html) with the following command:

```
python -m venv --system-site-packages pyEnv_name

source pyEnv_name/bin/activate

pip install --upgrade pip

pip install numpy networkx
```

## Validation of _GBA centrality_

All code for the validation of GBA centrality is in [Validation/](Validation/). For validation we used Python 3.12.