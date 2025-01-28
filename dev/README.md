### Part I: Leave-one-out validation

This is a comparison between the distribution of scores for left-out proteins with the distribution of scores for all proteins in the interactome. We compare the two distributions using a statistical test (Wilcoxon rank-sum) and calculate the p-value.

Example usage for calculations:

`python scripts/leave_one_out.py  -i input/Interactome_human.sif --causal_genes_file input/causalGenes.tsv --inUniProt /input/Uniprot_output.tsv --alpha 0.1 --d_max 10 --patho MMAF 1>leave_one_out/MMAF/alpha01_dmax10/scores_leave_out_out.tsv 2>leave_one_out/MMAF/alpha01_dmax10/log_leave_one_out.tsv`

### Part II: Tissue-enrichment validation

This is a comparison between the ratio of predicted causal proteins enriched in the tissue with the ratio of all proteins enriched in the tissue. We compare the two ratio using a statistical test (Fisher exact) and calculate the p-value to answer the question: "Are predicted causal genes significantly enriched in the tissue-specific genes?".

Expression Atlas was downloaded from Ensembl reference (v104) (https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv)

`mv E-MTAB-5214-query-results.tpms.tsv dev/`

Then add manually a column (after the ENSG column) "tissue ratio", which corresonds to the tissue-enrichment of each gene (i.e., expression_tissue / average_expression_tissue_all_genes).

### Part III: Robustness

This is a comparison between the lists of top 10% of highest-scoring genes generated with the new interactome (19/09/2024) vs the old interactome (15/01/2024). We hypothesize that the method is robust, therefore the two lists should be similar, even though the interactomes differ in the number of proteins and interactions. To determine the ovelap between the two lists, we compare then lists with a statistical test (Wilcoxon rank-sum) and calculate the p-value.

Use old interactome: input/Interactome_human_15012024.sif

`python scripts/GBA_centrality_PR.py -i input/Interactome_human_15012024.sif --causal_genes_file input/causalGenes_infertility.tsv --inUniProt input/Uniprot_output.tsv --alpha 0.5 --d_max 10 --patho MMAF 1>output/Interactome_human_15012024/MMAF/alpha05_d10/scores.tsv 2>output/Interactome_human_15012024/MMAF/alpha05_d10/log.txt`

Compare the overlap between the top 10% of highest-scoring genes from old interactome (15/01/2024) and new interactome (19/09/2024).
