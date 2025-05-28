############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024-2025
#
# This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################

import os
import sys
import logging
import pathlib

import argparse

import numpy
import networkx

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def calculate_scores(interactome, adjacency_matrices, causal_genes, alpha) -> dict:
    '''
    Calculate scores for every gene in the interactome based on the proximity to known causal genes.

    arguments:
    - interactome: type=networkx.Graph
    - adjacency_matrices: list of numpy arrays as returned by get_adjacency_matrices()
    - causal_genes: dict of causal genes with key=ENSG, value=1
    - alpha: attenuation coefficient (parameter set by user)

    returns:
    - scores: dict with key=ENSG, value=score
    '''
    # 1D numpy array for genes in the interactome: 1 if causal gene, 0 otherwise
    # size=len(nodes in interactome), ordered as in interactome.nodes()
    causal_genes_vec = numpy.zeros(len(interactome.nodes()), dtype=numpy.uint8)
    ni = 0
    for n in interactome.nodes():
        if n in causal_genes:
            causal_genes_vec[ni] = 1
        ni += 1

    # calculate scores
    scores_vec = numpy.zeros(len(causal_genes_vec))
    for d in range(0, len(adjacency_matrices)):
        A = adjacency_matrices[d]
        scores_vec += alpha ** d * A.dot(causal_genes_vec)

    scores = dict(zip(interactome.nodes(), scores_vec))  # map scores to genes

    return scores


def get_adjacency_matrices(interactome, d_max):
    '''
    Calculates powers of adjacency matrix up to power d_max (using non-normalized matrices).
    Then zeroes the diagonal and normalizes each matrix (row-wise).

    arguments:
    - interactome: type=networkx.Graph
    - d_max: max distance from a causal gene for it to contribute to a node's score (parameter set by user)

    returns:
    - adjacency_matrices: list of numpy arrays,
        array at index i (starting at i==0) is the processed A**i,
        rows and columns are ordered as in interactome.nodes()
    '''
    # list of processed matrices, should all be of the same type (eg numpy.float32)
    adjacency_matrices = []

    A = networkx.to_scipy_sparse_array(interactome, dtype=numpy.uint64, format='csc')

    # temporary variable for A**power
    res = numpy.identity(A.shape[0], dtype=numpy.uint64)
    adjacency_matrices.append(res.astype(numpy.float32))

    for power in range(1, d_max + 1):
        logger.debug(f"Calculating A**{power}")
        # @ - matrix multiplication, result of numpy @ CSC is numpy
        # res in CSR and A in CSC format should be fast, but it is slow,
        # also res and A in numpy format is very slow;
        # res in numpy format and A in CSC format is the fastest
        res = res @ A

        # zero the diagonal
        numpy.fill_diagonal(res, val=0)

        # row-wise normalization
        res_norm = res.astype(numpy.float32)
        row_sum = res_norm.sum(axis=1)  # row-wise axis=1, column-wise axis=0
        row_sum[row_sum == 0] = 1
        res_norm_T = res_norm.T
        res_norm_T /= row_sum

        adjacency_matrices.append(res_norm)

    logger.debug("Done building %i matrices", len(adjacency_matrices))
    return adjacency_matrices


def main(interactome_file, causal_genes_file, uniprot_file, alpha, d_max):

    logger.info("Parsing interactome")
    interactome = data_parser.parse_interactome(interactome_file)

    logger.info("Parsing gene-to-ENSG mapping")
    ENSG2gene, gene2ENSG, uniprot2ENSG = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, interactome)

    logger.info("Calculating powers of adjacency matrix")
    adjacency_matrices = get_adjacency_matrices(interactome, d_max)

    logger.info("Calculating scores")
    scores = calculate_scores(interactome, adjacency_matrices, causal_genes, alpha)

    logger.info("Printing scores")
    data_parser.scores_to_TSV(scores, ENSG2gene)

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
        prog=script_name,
        description="""
        GBA centrality is a network propagation algorithm for disease gene prioritization.
        The method assigns scores to genes that represent their likelihood of being causal for
        the phenotype/disease of interest. It takes into account the topology of the protein-protein
        interaction network (interactome) and prior knowledge about genes known to be associated with the disease.
        """
    )

    parser.add_argument('--interactome',
                        help='interactome SIF file with columns: ENSG1, "pp", ENSG2',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--causal',
                        help='TXT file (without a header) with 1 column: gene_name',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--uniprot',
                        help='parsed Uniprot file from Interactome/uniprot_parser.py',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--alpha',
                        help='attenuation coefficient (0 < alpha < 1)',
                        default=0.5,
                        type=float)
    parser.add_argument('--dmax',
                        help='propagation distance',
                        default=5,
                        type=int)

    args = parser.parse_args()

    try:
        main(interactome_file=args.interactome,
             causal_genes_file=args.causal,
             uniprot_file=args.uniprot,
             alpha=args.alpha,
             d_max=args.dmax)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
