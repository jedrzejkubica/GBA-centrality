### Part I: Leave-one-out validation

Example usage for calculations:

`python leave_one_out.py  -i ../input/Interactome_human.sif --causal_genes_file ../input/causalGenes.tsv --gene2ENSG_file ../input/canonicalGenes.tsv 1>leave_one_out/MMAF/alpha01_dmax10_alphanorm_10/scores_leave_out_out.tsv 2>leave_one_out/MMAF/alpha01_dmax10_alphanorm_10/log_leave_one_out.tsv --alpha 0.1 --d_max 10 --alpha_norm 10 --patho MMAF`

### Part II: Tissue-enrichment validation

### Part III: Contributions to c_i and N_i at various distances

Example usage for calculations: 

`python contributions_at_distances.py -i ../input/Interactome_human.sif --causal_genes_file ../input/causalGenes.tsv --gene2ENSG_file ../input/canonicalGenes.tsv --alpha 0.1 --d_max 10 --alpha_norm 10 --patho MMAF`

### Part IV: Robustness

Use old interactome: input/Interactome_human_15012024.sif
Output: ../output

`python scripts/GBA_centrality_PR.py -i input/Interactome_human_15012024.sif --causal_genes_file input/causalGenes_infertility.tsv --gene2ENSG_file input/canonicalGenes.tsv --alpha 0.5 --d_max 10 --patho MMAF  1> ~/workspace/output/Interactome_human_15012024/MMAF/alpha05_d10/scores.tsv 2> ~/workspace/output/Interactome_human_15012024/MMAF/alpha05_d10/log.txt`

Compare top 10% of highest-scoring genes from old interactome (15_01_2024) and new interactome (19_09_2024)
