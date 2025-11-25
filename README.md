> Laboratory TIMC / BCM Grenoble

# GBA (Guilt-by-association) centrality

GBA centrality is a network propagation algorithm for disease gene prioritization. The method assigns scores to genes that represent their likelihood of being causal for the phenotype/disease of interest. It takes into account the topology of the protein-protein interaction network (interactome) and prior knowledge about genes known to be associated with the disease.

## How to use GBA centrality

We assume here that this repository is cloned into `~/Software/`, input data is downloaded and processed in `~/GBA-input/`, and results will be produced in `~/GBA-output/`. Change these names to your taste and adapt all commands below accordingly. This repository requires [GBA-centrality-C](https://github.com/jedrzejkubica/GBA-centrality-C), because `GBA_centrality.py` uses a GBA-centrality-C shared object (.so file) for heavy-lifting calculations. Create these folders and set up GBA-centrality-C with the following commands:

```
mkdir ~/Software/ ~/GBA-input/ ~/GBA-output/

cd ~/Software/
git clone https://github.com/jedrzejkubica/GBA-centrality.git
cd GBA-centrality
git submodule init
git submodule update
cd GBA-centrality-C
make
```

As input, GBA centrality takes an interactome SIF file, a parsed Uniprot file and a file with known disease-associated genes. GBA centrality allows the user to set attenuation coefficient `--alpha`  (0 < alpha < 1 ; default = 0.5).

Example usage:
```
python ~/Software/GBA-centrality/GBA_centrality.py \
  --interactome ~/GBA-input/interactome_human.sif \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  --causal ~/GBA-input/causal_genes.txt \
  --alpha 0.5 \
  1> ~/GBA-output/scores.tsv \
  2> ~/GBA-output/log.txt
```


## How to prepare input data

### Uniprot file

This file is used for mapping between gene names, gene ENSGs and protein Uniprot IDs.

Download and parse Uniprot (file size ~600Mb):

```
cd ~/GBA-input
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
gunzip -c uniprot_sprot.dat.gz | python ~/Software/GBA-centrality/Interactome/uniprot_parser.py > uniprot_parsed.tsv
```

### Interactome SIF file

Below we provide commands for constructing a human interactome (undirected and unweighted) using protein-protein interaction data from [BioGRID](https://thebiogrid.org/), [IntAct](https://www.ebi.ac.uk/intact/home) and [Reactome](https://reactome.org/download-data).

**Step 1. Download and extract human protein-protein interaction data**

BioGRID (file size ~170Mb)

```
cd ~/GBA-input
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
unzip BIOGRID-ORGANISM-LATEST.mitab.zip BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt
```

IntAct (file size ~800Mb)

```
cd ~/GBA-input
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
```

Reactome (file size ~170Mb)

```
cd ~/GBA-input
wget https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.psi-mitab.txt
```

**Step 2. Parse protein-protein interaction data**

Parse BioGRID

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  > ~/GBA-input/interactions_Biogrid.tsv
```

Parse IntAct

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/intact.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  > ~/GBA-input/interactions_Intact.tsv
```

Parse Reactome

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/reactome.homo_sapiens.interactions.psi-mitab.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  > ~/GBA-input/interactions_Reactome.tsv
```

**Step 3. Build a high-confidence human interactome**

Here the protein-protein interactions data is filtered on "Interaction Detection Method" and "Interaction Type": each interaction must be confirmed by at least 1 experiment using a binary interaction detection method. We also remove the "artificial protein hubs" (i.e., top 1% of proteins with the highest degrees). The data is merged into one SIF file (for details about the SIF format see: https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format).

```
python ~/Software/GBA-centrality/Interactome/build_interactome.py \
  --interactions ~/GBA-input/interactions_Biogrid.tsv ~/GBA-input/interactions_Intact.tsv ~/GBA-input/interactions_Reactome.tsv \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  > ~/GBA-input/interactome_human.sif
```

The interactome SIF file has 3 tab-separated columns: ENSG1 "pp" ENSG2. It can be modified to be directed and/or weighted: for a weighted graph replace "pp" with weights (0 < weight < 1) and use `--weighted` parameter; for a directed graph use `--directed` parameter. Then run GBA centrality as follows:

```
python ~/Software/GBA-centrality/GBA_centrality.py --weighted --directed [...]
```


### File with known disease-associated genes

Create a file `causal_genes.txt` (without a header) with one known causal gene name per line.

> [!NOTE]
> GBA centrality maps causal gene names to ENSGs using the parsed Uniprot file: As gene names, it requires the HGNC nomenclature (HUGO Gene Nomenclature Committee, https://www.genenames.org).


## Validation of GBA centrality

All code for the validation of GBA centrality can be found here: [GBA-centrality-validation](https://github.com/jedrzejkubica/GBA-centrality-validation).
