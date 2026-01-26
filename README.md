> Laboratory TIMC / BCM Grenoble

# Guilt-by-association centrality

GBA centrality is a new network propagation algorithm based on non-backtracking walks and in-degree normalization. The method assigns scores to nodes in the network that represent their likelihood of being in proximity to a given list of nodes of interest.

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

As input, GBA centrality takes a network SIF file (for details about the SIF format see: https://cytoscape.org/manual/Cytoscape2_5Manual.html#SIF%20Format) and a file with nodes of interest (seeds). GBA centrality allows the user to set attenuation coefficient `--alpha`  (0 < alpha < 1 ; default = 0.5).

Example usage:
```
python ~/Software/GBA-centrality/GBA_centrality.py \
  --network ~/GBA-input/interactome_human.sif \
  --seeds ~/GBA-input/causal_genes.txt \
  --alpha 0.5 \
  1> ~/GBA-output/scores.tsv \
  2> ~/GBA-output/log.txt
```



## Validation of GBA centrality

All code for the validation of GBA centrality can be found here: [GBA-centrality-validation](https://github.com/jedrzejkubica/GBA-centrality-validation).
