## Contents:

- [validation_GBA_centrality.ipynb](validation_GBA_centrality.ipynb)

    For a complete description of the validation of GBA centrality see: (ref)

    > [!NOTE]
    > For tissue enrichment validation download the Expression Atlas from Ensembl reference (v104)
    > (https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv).
    > Then add manually a column (after the ENSG column) called "tissue_enrichment", which corresonds to the tissue enrichment
    > of each gene (i.e., divide tissue expression of each gene by the average expression of that gene in all tissues).

- [leave_one_out_scores.py](leave_one_out_scores.py)

    This script calculates the scores for left-out genes.

    example usage:
    ```
    python leave_one_out_scores.py \
        -i input/interactome_human.sif \
        --causal_genes_file input/causal_genes_infertility.tsv \
        --uniprot_file input/uniprot_output.tsv \
        --patho MMAF \
        --alpha 0.5 \
        --d_max 10 \
        1> output/scores_leave_one_out.tsv \
        2> output/log_leave_one_out_scores.txt
    ```

- [leave_one_out_ranks.py](leave_one_out_ranks.py)

    This script is similar to leave_one_out_scores.py but prints to STDOUT ranks for left-out genes instead of scores.