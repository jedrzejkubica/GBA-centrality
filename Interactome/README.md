> Laboratory TIMC / BCM Grenoble

## How to build a human interactome




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

**Step 3. Build a human interactome**

Here the interactome data is filtered on "Interaction Detection Method" and "Interaction Type": each interaction must be confirmed by at least 1 experiment using a binary interaction detection method. The data is merged into one SIF file.

```
python ~/Software/GBA-centrality/Interactome/build_interactome.py \
  --interactions ~/GBA-input/interactions_Biogrid.tsv ~/GBA-input/interactions_Intact.tsv ~/GBA-input/interactions_Reactome.tsv \
  > ~/GBA-input/interactome_human.sif
```

The interactome SIF file has 3 tab-separated columns: protein1 "pp" protein2.


### File with seeds

known disease-associated genes

Create a file `causal_genes.txt` (without a header) with one known causal gene name per line.

> [!NOTE]
> GBA centrality maps causal gene names to ENSGs using the parsed Uniprot file: As gene names, it requires the HGNC nomenclature (HUGO Gene Nomenclature Committee, https://www.genenames.org).


`causal_gene_parser.py` will take `causal_genes.txt` and produce a list of Primary AC for proteins corresponding to the causal genes in `seeds.txt`.

```
python ~/Software/GBA-centrality/Interactome/causal_genes_parser.py \
  --uniprot ~/GBA-input/uniprot_parsed.tsv \
  --causal ~/GBA-input/causal_genes.txt
  > ~/GBA-input/causal_proteins.txt
```