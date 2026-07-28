> Laboratory TIMC / BCM Grenoble


# Backtrack-Free Walk (BFWalk)

BFWalk is a new network propagation algorithm based on non-backtracking walks and in-degree normalization. The algorithm assigns scores to nodes in the network that represent their proximity to a given list of nodes of interest (i.e. seeds).


## Install BFWalk

This repository requires [BFWalk-C](https://github.com/jedrzejkubica/BFWalk-C), because `BFWalk.py` uses a BFWalk-C shared object (.so file) for heavy-lifting calculations. Set up BFWalk with the following commands:

```
git clone https://github.com/jedrzejkubica/BFWalk.git
cd BFWalk
git submodule init
git submodule update
cd BFWalk-C
make
cd ..
```


## Use BFWalk

For details see:
```
python BFWalk.py --help
```

As input, BFWalk requires:

- `--network`: a text file with one interaction per line, in the following format: `node1 weight/interaction_type node2` (3 tab-separated columns). It is a format similar to SIF but allows for weighted networks ([SIF format documentation](https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format))
- `--seeds`: a text file with one seed per line

For output, BFWalk prints to stdout the scores in TSV format `node score` (2 tab-separated columns).

Networks are undirected and unweighted by default, but BFWalk can also use a directed and/or weighted network:
- use `--directed` if your network is directed, edges are then seen as node1->node2 (i.e. the source is in the first column and the destination in the third column);
- use `--weighted` if your network is weighted, the second column of the network file must then contain the weight of each interaction (0 < weight <= 1).

If needed, BFWalk allows the user to set the attenuation coefficient `--alpha`  (0 < alpha < 1), although the default = 0.5 should be fine for most use cases.


## Examples

### Weighted network

This example uses a simple "diamond" network with 4 nodes and 4 weighted edges: A, B, C, D. Here A is the seed.

```
python BFWalk.py \
  --network Examples/network_weighted.sif \
  --seeds Examples/seeds.txt \
  --weighted \
  1> scores.tsv \
  2> log.txt
```


### Directed network

This example uses a simple network with 3 nodes and 2 directed edges: C -> A -> B. Here A is the seed.

```
python BFWalk.py \
  --network Examples/network_directed.sif \
  --seeds Examples/seeds.txt \
  --directed \
  1> scores.tsv \
  2> log.txt
```


### Human interactome

We provide additional instructions for interactome-based disease gene prioritization. The instructions and scripts to build a human interactome and prepare a seeds file with known causal genes can be found in the [Interactome/](Interactome/) subdirectory.


## Validation of BFWalk

A manuscript describing BFWalk has been submitted. The code to perform the analyses and generate the figures presented in this manuscript are available on GitHub: [BFWalk-validation](https://github.com/jedrzejkubica/BFWalk-validation).


## Note

Formerly "GBA-centrality" for anyone arriving via old citations.