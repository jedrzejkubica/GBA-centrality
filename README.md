> Laboratory TIMC / BCM Grenoble

# Guilt-by-association centrality

GBA centrality is a new network propagation algorithm based on non-backtracking walks and in-degree normalization. The method assigns scores to nodes in the network that represent their proximity to a given list of nodes of interest (i.e. seeds).

## Install GBA centrality

This repository requires [GBA-centrality-C](https://github.com/jedrzejkubica/GBA-centrality-C), because `GBA_centrality.py` uses a GBA-centrality-C shared object (.so file) for heavy-lifting calculations. Set up GBA-centrality-C with the following commands:

```
git clone https://github.com/jedrzejkubica/GBA-centrality.git
cd GBA-centrality
git submodule init
git submodule update
cd GBA-centrality-C
make
cd ..
```


## Use GBA centrality

As input, GBA centrality takes:

- network: a text file with one interaction per line, in the following format: source weight/interaction_type destination (3 tab-separated columns). It is a format similar to SIF but allows for weighted networks (for details about the SIF format see: https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format)
- seeds: a text file with one seed per line

GBA centrality can also use a directed and/or weighted network: for a weighted network put weights (0 < weight <= 1) in the second column of the network file and use `--weighted` parameter; for a directed network use `--directed` parameter. Then run GBA centrality as follows:

```
python GBA_centrality.py --weighted --directed [...]
```

If needed, GBA centrality allows the user to set the attenuation coefficient `--alpha`  (0 < alpha < 1), although the default = 0.5 should be fine for most use cases.

Example usage:
```
python GBA_centrality.py \
  --network ~/GBA-input/network.sif \
  --seeds ~/GBA-input/seeds.txt \
  --alpha 0.5 \
  1> ~/GBA-output/scores.tsv \
  2> ~/GBA-output/log.txt
```

## Examples

### Weighted network

This example uses a simple "diamond" network with 4 nodes and 4 weighted edges: A, B, C, D. Here A is the seed.

Example usage:
```
python GBA_centrality.py \
  --network Example/network_weighted.sif \
  --seeds Example/seeds.txt \
  --weighted
```


### Directed network

This example uses a simple network with 3 nodes and 2 directed edges: C -> A -> B. Here A is the seed.

Example usage:
```
python GBA_centrality.py \
  --network Example/network_directed.sif \
  --seeds Example/seeds.txt \
  --directed
```


### Human interactome

The scripts used to build a human interactome and seeds file with causal genes/proteins can be found here: [Interactome/](Interactome/). For detailed instructions see here [Interactome/README.md](Interactome/README.md).

All code to perform the analyses and generate the figures are also available on GitHub: [GBA-centrality-validation](https://github.com/jedrzejkubica/GBA-centrality-validation).
