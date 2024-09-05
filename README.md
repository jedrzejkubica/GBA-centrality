> PhD project at TIMC / MAGe Grenoble (September 2024 - )

# Interactome-TIMC

This repository contains scripts for _GBA (Guilt-by-association) centrality_. 

The method assigns a score for every protein in the interactome that represents the likelihood of its functional association with the phenotype of interest.


## Example usage of _GBA centrality_

As input, it requires an interactome, canonical genes, as well as known causal genes for the phenotype of interest.

Example usage for infetility MMAF phenotype and parameters alpha=0.1, d_max=10, alpha_norm=10:

```
python scripts/GBA_centrality.py \
  -i input/Interactome_human.sif \
  --causal_genes_file input/causalGenes.tsv \
  --gene2ENSG_file input/canonicalGenes.tsv \
  --alpha 0.1 \
  --d_max 10 \
  --alpha_norm 10 \
  --patho MMAF \
  1>output/scores.tsv \
  2>output/log.txt
```


## Prepare input data

### Step 1. Download interactome data

Download Uniprot data

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

```
gunzip -c uniprot_sprot.dat.gz
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
python scripts/Uniprot_parser.py > Uniprot_output.tsv
```

Parse BioGRID

```
python scripts/Interaction_parser.py --inInteraction BIOGRID-ORGANISM-Homo_sapiens*.mitab.txt --inUniProt Uniprot_output.tsv > Exp_Biogrid.tsv
```

Parse IntAct

```
python scripts/Interaction_parser.py --inInteraction intact.txt --inUniProt Uniprot_output.tsv > Exp_Intact.tsv
```

Parse Reactome

```
python scripts/Interaction_parser.py --inInteraction reactome.homo_sapiens.interactions.psi-mitab.txt --inUniProt Uniprot_output.tsv > Exp_Reactome.tsv
```


### Step 3. Build the interactome

First create a TSV file `canonicalGenes.tsv` (with a header "GENE\tENSG") with 2 columns: gene_name, ENSG.

> [!NOTE]
> We provide an example file with canonical genes in [input/canonicalGenes.tsv](input/canonicalGenes.tsv).

Then build the interactome:

```
python scripts/Build_Interactome.py \
  --inExpFile Exp_Biogrid.tsv Exp_Intact.tsv Exp_Reactome.tsv \
  --inUniProt Uniprot_output.tsv \
  --inCanonicalFile canonicalGenes.tsv > Interactome_human.tsv
```
> [!NOTE]
> We provide an example human interactome in [input/Interactome_human.sif](input/Interactome_human.sif) (with 14,393 proteins and 89,319 interactions).


### Step 4. Prepare a list of causal genes

Create a TSV file `causalGenes.tsv` (without a header) with 2 columns: gene_name, pathologyID.

> [!NOTE]
> We provide an example file with causal genes for infertility phenotypes (e.g., MMAF) in [input/causalGenes.tsv](input/causalGenes.tsv).


## Validation of _GBA centrality_

We put scripts and jupyter notebooks for validation in [dev/](dev/).


### Python environment

We used Python 3.9 and provided dependencies in [requirements.txt](requirements.txt) (venv; https://docs.python.org/3/library/venv.html).
