> Laboratory TIMC / BCM Grenoble

# GBA (Guilt-by-association) centrality

GBA centrality is a network propagation algorithm for disease gene prioritization. The method assigns scores to genes that represent their likelihood of being causal for the phenotype/disease of interest. It takes into account the topology of the protein-protein interaction network (interactome) and prior knowledge about genes known to be associated with the disease.

## 🚀 How to use _GBA centrality_

We assume here that code is cloned in ~/Software/ , input data is downloaded and processed in ~/GBA-input/ , and results are produced in ~/GNA-output/ . Change these names to your taste and adapt all commands below accordingly. Create them with the following command:
```
mkdir ~/Software/ ~/GBA-input/ ~/GNA-output/
```

As input, GBA centrality takes an interactome SIF file, a parsed Uniprot DAT file and a TSV file with known disease-associated genes.

Example usage for an infertility phenotype (MMAF: multiple morphological abnormalities of the sperm flagella):
```
python ~/Software/GBA-Centrality/GBA_centrality.py \
  -i ~/GBA-input/interactome_human.sif \
  --uniprot_file ~/GBA-input/uniprot_parsed.tsv \
  --causal_genes_file ~/GBA-input/causal_genes_infertility.tsv \
  --patho MMAF \
  1> ~/GBA-output/scores.tsv \
  2> ~/GBA-output/log.txt
```

### Parameters

GBA centrality allows the user to set two parameters:
- alpha (default = 0.5) - attenuation coefficient; 0 < alpha < 1
- d_max (default = 5) - propagation distance; max. distance for a causal gene to contibute to a gene's score

Example:

```
python GBA_centrality.py --alpha 0.5 --d_max 5 [...]
```

## How to prepare input data

### Uniprot DAT file

Download Uniprot data

```
cd ~/GBA-input
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

Parse uniprot_sprot.dat.gz

```
gunzip -c uniprot_sprot.dat.gz | python Interactome/uniprot_parser.py > uniprot_parsed.tsv
```

### Interactome SIF file

**Step 1. Download and extract Human protein-protein interaction data**

[BioGRID](https://thebiogrid.org/)

```
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
unzip BIOGRID-ORGANISM-LATEST.mitab.zip BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt
```

[IntAct](https://www.ebi.ac.uk/intact/home)

```
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
```

[Reactome](https://reactome.org/download-data)

```
wget https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.psi-mitab.txt
```

**Step 2. Parse protein-protein interaction data**

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

**Step 3. Build a high-confidence human interactome**

```
python Interactome/build_interactome.py \
  --interactions_parsed_files interactions_Biogrid.tsv interactions_Intact.tsv interactions_Reactome.tsv \
  --uniprot_file uniprot_parsed.tsv > interactome_human.sif
```

> [!NOTE]  
> The build_interactome.py script maps protein Uniprot IDs to gene ENSG IDs.

### TSV file with known disease-associated genes

Create a tab-separated file `causalGenes.tsv` (without a header) with 2 columns: gene_name, pathology

> [!NOTE]
> GBA centrality maps disease-assocaited gene names to ENSG IDs using the Uniprot DAT file.
> 
> As gene names, GBA centrality requires the HGNC nomenclature (HUGO Gene Nomenclature Committee, https://www.genenames.org).

### Python environment

GBA centrality is written in Python :snake: and requires the following dependencies: [NumPy](https://numpy.org/) and [NetworkX](https://networkx.org/).

We recommend installing them via [Python venv](https://docs.python.org/3/library/venv.html) with the following command:

```
python -m venv --system-site-packages ~/pyEnv_GBA-centrality
source ~/pyEnv_GBA-centrality/bin/activate
pip install --upgrade pip
pip install numpy networkx
```

You can then run GBA-centrality with:
```
source ~/pyEnv_GBA-centrality/bin/activate
python GBA_centrality.py [...]
```

## Validation of _GBA centrality_

All code for the validation of GBA centrality is in [GBA-centrality-validation](https://github.com/jedrzejkubica/GBA-centrality-validation). For validation we used Python 3.12.
