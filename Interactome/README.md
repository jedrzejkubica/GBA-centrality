> Laboratory TIMC / BCM Grenoble

## Build a human interactome

Below we provide instructions for interactome-based disease gene prioritization. We provide scripts and commands to: download and parse protein-protein interaction data, build a human interactome, and create a file with seeds (i.e. disease-associated genes/proteins).

We assume that the `GBA-centrality/` repository is in `~/Software/`, and that data will be downloaded into `~/GBA-input`. If needed, adapt the commands below to your taste.

```
mkdir ~/GBA-input
cd ~/GBA-input
```

#### Uniprot file

This file will be used for for mapping between Uniprot ACs, tax IDs and gene names.

Download and parse Uniprot (file size ~660Mb):

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
gunzip -c uniprot_sprot.dat.gz | python ~/Software/GBA-centrality/Interactome/uniprot_parser.py > uniprot_parsed.tsv
```


### Interactome SIF file

Build a human interactome (undirected and unweighted) using protein-protein interaction (PPI) data from [BioGRID](https://thebiogrid.org/) and [IntAct](https://www.ebi.ac.uk/intact/home).

**Step 1. Download and extract human PPI data**

BioGRID (zipped file size ~176Mb)

```
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
unzip BIOGRID-ORGANISM-LATEST.mitab.zip BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt
```

IntAct (zipped file size ~1.3Gb)

```
wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
```


**Step 2. Parse PPI data**

The PPI data will be filtered on "Interaction Detection Method", ie. ignore some "bad" detection methods (MI:0254 (genetic interference) and MI:0686 (unspecified method)).

Parse BioGRID

```
python ~/Software/GBA-centrality/Interactome/interaction_parser.py \
  --interactions ~/GBA-input/BIOGRID-ORGANISM-Homo_sapiens\*.mitab.txt \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  1> ~/GBA-input/interactions_Biogrid.tsv \
  2> ~/GBA-input/interactions_Biogrid.log
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

The interactome will be saved in a file format similar to SIF (https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format)
with 3 tab-separated columns: protein1 "pp" protein2.

```
python ~/Software/GBA-centrality/Interactome/build_interactome.py \
  --interactions ~/GBA-input/interactions_Biogrid.tsv ~/GBA-input/interactions_Intact.tsv \
  > ~/GBA-input/interactome_human.sif
```

If needed, `build_interactome.py` allows the user to set the min number of evidences `--n_evidence`  (default=2) and the min number of direct interactions `--n_direct` (default=1).


### Seeds file

Seeds must correspond to nodes in the interactome. Because the interactome contains Uniprot ACs, causal genes provided as seeds must be mapped to Uniprot Primary ACs beforehand. We provide the script `causal_genes_parser.py`, which maps causal gene names to UniProt Primary ACs using the `uniprot_parsed.tsv` file.

Create a file `causal_genes.txt` (without a header) with one causal gene name per line.

```
python ~/Software/GBA-centrality/Interactome/causal_genes_parser.py \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  --causal ~/GBA-input/causal_genes.txt \
  > ~/GBA-input/causal_proteins.txt
```

Causal proteins will be saved in `causal_proteins.txt`, one UniProt Primary AC per line. It can happend that some genes are mapped to more than one protein when there are multiple, genuinely distinct protein products (in that case all proteins will be saved as causal).
