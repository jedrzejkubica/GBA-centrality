> Laboratory TIMC / BCM Grenoble

## Build a human interactome

Below we provide instructions for interactome-based disease gene prioritization. We provide scripts and commands to: download and parse protein-protein interaction data, build a human interactome, and create a file with seeds (i.e. disease-associated genes/proteins).

We assume that the `GBA-centrality/` repository is in `~/Software/`, and that data will be downloaded into `~/GBA-input`. If needed, adapt the commands below to your taste.

```
mkdir ~/GBA-input
cd ~/GBA-input
```

#### Uniprot file

This file will be used for for mapping between gene names, gene ENSGs and protein Uniprot IDs.

Download and parse Uniprot (file size ~660Mb):

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
gunzip -c uniprot_sprot.dat.gz | python ~/Software/GBA-centrality/Interactome/uniprot_parser.py > uniprot_parsed.tsv
```


### Interactome SIF file

Build a human interactome (undirected and unweighted) using protein-protein interaction (PPI) data from [BioGRID](https://thebiogrid.org/), [IntAct](https://www.ebi.ac.uk/intact/home) and [Reactome](https://reactome.org/download-data).

**Step 1. Download and extract human PPI data**

BioGRID Multi-Validated (MV) Datasets (file size ~35Mb)

```
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.mitab.zip
unzip BIOGRID-MV-Physical-LATEST.mitab.zip
```

IntAct (file size ~1.3Gb)

```
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
```


**Step 2. Parse PPI data**

Parse BioGRID

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/BIOGRID-MV-Physical-*.mitab.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  1> ~/GBA-input/interactions_Biogrid-MV.tsv \
  2> ~/GBA-input/interactions_Biogrid-MV.log
```

Parse IntAct

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/intact.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  1> ~/GBA-input/interactions_Intact.tsv \
  2> ~/GBA-input/interactions_Intact.log
```


**Step 3. Build a human interactome**

The interactome data will be filtered on "Interaction Detection Method" and "Interaction Type": each interaction must be confirmed by at least 1 binary detection method. Then the interactome will be saved in a file format similar to SIF (https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format).

```
python ~/Software/GBA-centrality/Interactome/build_interactome.py \
  --interactions ~/GBA-input/interactions_Biogrid-MV.tsv ~/GBA-input/interactions_Intact.tsv \
  > ~/GBA-input/interactome_human.sif
```

If needed, `build_interactome.py` allows the user to set the min number of evidences `--n_evidence`  (default=2) and the min number of direct interactions `--n_direct` (default=1).

The interactome file has 3 tab-separated columns: protein1 "pp" protein2.


### Seeds file

Seeds must correspond to nodes in the interactome. Because the interactome contains protein IDs, causal genes provided as seeds must be mapped to protein IDs beforehand. We provide the script `causal_genes_parser.py`, which maps causal gene names to UniProt accession numbers (UniProt Primary ACs) using the parsed Uniprot file.

> [!NOTE]
> It requires [HUGO Gene Nomenclature Committee](https://www.genenames.org) gene names.

Create a file `causal_genes.txt` (without a header) with one causal gene name per line. Then map gene names in `causal_genes.txt` into protein IDs and save them in `causal_proteins.txt`.


```
python ~/Software/GBA-centrality/Interactome/causal_genes_parser.py \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  --causal ~/GBA-input/causal_genes.txt
  > ~/GBA-input/causal_proteins.txt
```
