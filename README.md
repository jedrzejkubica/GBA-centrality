> Laboratory TIMC / BCM Grenoble

# GBA (Guilt-by-association) centrality

GBA centrality is a network propagation algorithm for disease gene prioritization. The method assigns scores to genes that represent their likelihood of being causal for the phenotype/disease of interest. It takes into account the topology of the protein-protein interaction network (interactome) and prior knowledge about genes known to be associated with the disease.

## 🚀 How to use GBA centrality

We assume here that this repository is cloned into `~/Software/`, input data is downloaded and processed in `~/GBA-input/`, and results will be produced in `~/GBA-output/`. Change these names to your taste and adapt all commands below accordingly. Create these folders with the following command:

```
mkdir ~/Software/ ~/GBA-input/ ~/GBA-output/
```

```
git clone https://github.com/jedrzejkubica/GBA-centrality.git ~/Software/
```

This repository uses [GBA-centrality-C](https://github.com/jedrzejkubica/GBA-centrality-C). GBA_centrality.py uses the GBA-centrality-C shared object (.so file) created using the following commands:

```
cd GBA-centrality/
git config submodule.GBA-centrality-C.url git@github.com:jedrzejkubica/GBA-centrality-C.git
git submodule init
git submodule update
(for development) git checkout submodule-in-C ; git submodule update
cd GBA-centrality-C
make
```

As input, GBA centrality takes an interactome SIF file, a parsed Uniprot file and a TXT file with known disease-associated genes.

Example usage:
```
python ~/Software/GBA-centrality/GBA_centrality.py \
  --interactome ~/GBA-input/interactome_human.sif \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  --causal ~/GBA-input/causal_genes.txt \
  1> ~/GBA-output/scores.tsv \
  2> ~/GBA-output/log.txt
```

### Parameters

GBA centrality allows the user to set alpha - attenuation coefficient (0 < alpha < 1 ; default = 0.5).

Example:
```
python ~/Software/GBA-centrality/GBA_centrality.py --alpha 0.5 [...]
```

## How to prepare input data

### Uniprot file

This file will be used for mapping between Uniprot protein IDs to gene ENSG IDs and gene names to gene ENSG IDs.

Download and parse Uniprot dataset (.gz file size ~651M):

```
cd ~/GBA-input
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
gunzip -c uniprot_sprot.dat.gz | python ~/Software/GBA-centrality/Interactome/uniprot_parser.py > uniprot_parsed.tsv
```

### Interactome SIF file

Below we provide commands for constructing a human interactome (undirected and unweighted) using protein-protein interaction datasets from [BioGRID](https://thebiogrid.org/), [IntAct](https://www.ebi.ac.uk/intact/home) and [Reactome](https://reactome.org/download-data). The interactome file has 3 tab-separated columns: ENSG1 "pp" ENSG2. It can be modified to create a weighted graph by replacing "pp" with weights (0 < weight < 1). If so, use `--weighted` parameter when running GBA centrality.

**Step 1. Download and extract human protein-protein interaction data**

BioGRID (.zip file size ~171M)

```
cd ~/GBA-input
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
unzip BIOGRID-ORGANISM-LATEST.mitab.zip BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt
```

IntAct (.zip file size ~807M)

```
cd ~/GBA-input
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
```

Reactome (file size ~176M)

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

Here the protein-protein interactions data is filtered as follows:
- remove complex expansion
- remove interaction detection methods: genetic interference (MI:0254) and unspecified method (MI:0686)
- keep interaction types: direct interaction (MI:0407) and physical association (MI:0915)
- keep interactions if at least one experiment has been proven by a binary interaction method
- remove interactions proven by Affintity Chromatography Technology (MI:0004)
- remove self-loops
- remove the top 1% of proteins with the highest degrees (number of interactions)

Then the data is merged into one SIF file. Further details about the SIF format can be found here: https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format

```
python ~/Software/GBA-centrality/Interactome/build_interactome.py \
  --interactions ~/GBA-input/interactions_Biogrid.tsv ~/GBA-input/interactions_Intact.tsv ~/GBA-input/interactions_Reactome.tsv \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  > ~/GBA-input/interactome_human.sif
```

> [!NOTE]
> The build_interactome.py script maps protein Uniprot IDs to gene ENSG IDs using the parsed Uniprot file.

### File with known disease-associated genes

Create a file `causal_genes.txt` (without a header) with one known causal gene name per line.

> [!NOTE]
> GBA centrality maps causal gene names to ENSG IDs using the parsed Uniprot file.
> 
> As gene names, GBA centrality requires the HGNC nomenclature (HUGO Gene Nomenclature Committee, https://www.genenames.org).

### Python environment

GBA centrality is written in Python :snake:

We recommend setting up a [Python venv](https://docs.python.org/3/library/venv.html) with the following command:

```
python -m venv --system-site-packages ~/pyEnv_GBA-centrality
source ~/pyEnv_GBA-centrality/bin/activate
pip install --upgrade pip
```

You can then run GBA-centrality with:
```
source ~/pyEnv_GBA-centrality/bin/activate
python ~/Software/GBA-centrality/GBA_centrality.py [...]
```

## Validation of GBA centrality

All code for the validation of GBA centrality is in [GBA-centrality-validation](https://github.com/jedrzejkubica/GBA-centrality-validation). For validation we used Python 3.12.
