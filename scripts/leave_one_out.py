import logging
import os
import sys

import pathlib

import argparse

import GBA_centrality_PR
import data_parser


def leave_one_out(interactome, adjacency_matrices, causal_genes, alpha):
    '''
    arguments:
    - interactome: type=networkx.Graph
    - adjacency_matrices: list of scipy sparse arrays as returned by GBA_centrality.get_adjacency_matrices()
    - causal_genes: dict of causal genes with key=ENSG, value=1

    returns:
    - ranks_left_out: dict with key=causal gene (ENSG), value=rank of this gene
      when it is left out
    '''
    # initialize dict to store left-out scores
    ranks_left_out = {}

    for left_out in list(causal_genes.keys()):
        logger.info("Leaving out %s", left_out)
        del causal_genes[left_out]
        scores = GBA_centrality_PR.calculate_scores(interactome, adjacency_matrices, causal_genes, alpha)

        # save the left-out rank
        results_sorted = sorted(scores.keys(), key=lambda item: scores[item], reverse=True)
        rank_left_out = results_sorted.index(left_out) + 1 # ranks start at 1

        ranks_left_out[left_out] = rank_left_out
        causal_genes[left_out] = 1

    return ranks_left_out


def main(interactome_file, causal_genes_file, Uniprot_file, patho, alpha, d_max):

    logger.info("Parsing interactome")
    interactome = data_parser.parse_interactome(interactome_file)

    logger.info("Parsing gene-to-ENSG mapping")
    ENSG2gene, gene2ENSG, Uniprot2ENSG = data_parser.parse_Uniprot(Uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, interactome, patho)

    logger.info("Calculating powers of adjacency matrix")
    adjacency_matrices = GBA_centrality_PR.get_adjacency_matrices(interactome, d_max)

    logger.info("Calculating leave-one-out ranks")
    ranks = leave_one_out(interactome, adjacency_matrices, causal_genes, alpha)

    logger.info("Printing leave-one-out ranks")
    data_parser.scores_to_TSV(ranks, ENSG2gene)

    logger.info("Done!")


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog="GBA_centrality.py",
        description="Calculate leave-one-out for GBA centrality."
    )

    parser.add_argument('-i', '--interactome_file', type=pathlib.Path, required=True)
    parser.add_argument('--causal_genes_file', type=pathlib.Path, required=True)
    parser.add_argument('--Uniprot_file', type=pathlib.Path, required=True)
    parser.add_argument('--patho', default='MMAF', type=str)
    parser.add_argument('--alpha', default=0.5, type=float)
    parser.add_argument('--d_max', default=5, type=int)

    args = parser.parse_args()

    try:
        main(interactome_file=args.interactome_file,
             causal_genes_file=args.causal_genes_file,
             Uniprot_file=args.Uniprot_file,
             patho=args.patho,
             alpha=args.alpha,
             d_max=args.d_max)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
