############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024
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

import logging
import networkx
import numpy
import os
import sys

import pathlib

import argparse

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def calculate_scores(interactome, adjacency_matrices, causal_genes, alpha) -> dict:
    '''
    Calculates scores for every gene in the interactome based on the proximity to causal genes.
    The algorithm:
    1. Calculate all A**k up to k == d_max
    2. Normalize all A**k row-wise
    2. Zero the diagonals of all A**k
    3. Multiply A**k and c (where 1: causal, 0: non-causal)

    arguments:
    - interactome: type=networkx.Graph
    - adjacency_matrices: list of numpy arrays as returned by get_adjacency_matrices()
    - causal_genes: dict of causal genes with key=ENSG, value=1
    - alpha: attenuation parameter

    returns:
    - scores: dict with key=ENSG, value=score
    '''
    # 1D numpy array for genes in the interactome: 1 if causal gene, 0 otherwise,
    # size=len(nodes in interactome), ordered as in interactome.nodes()
    causal_genes_vec = numpy.zeros(len(interactome.nodes()), dtype=numpy.uint8)
    ni = 0
    for n in interactome.nodes():
        if n in causal_genes:
            causal_genes_vec[ni] = 1
        ni += 1

    scores_vec = numpy.zeros(len(causal_genes_vec))

    # calculate scores
    for d in range(0, len(adjacency_matrices)):
        A = adjacency_matrices[d]
        scores_vec += alpha ** d * A.dot(causal_genes_vec)

    # map ENSGs to scores
    scores = dict(zip(interactome.nodes(), scores_vec))

    return scores


def get_adjacency_matrices(interactome, d_max=5):
    '''
    Calculates powers of adjacency matrix up to power d_max (using non-normalized matrices).
    Then zeroes the diagonal and normalizes the matrix (row-wise).

    arguments:
    - interactome: type=networkx.Graph
    - d_max: int

    returns:
    - adjacency_matrices: list of numpy arrays, array at index i (starting at i==0)
      is the processed A**i, rows and columns are ordered as in interactome.nodes()
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
        # also res and A in numpy format is very slow,
        # but res in numpy format and A in CSC format is the fastest
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


def main(interactome_file, causal_genes_file, Uniprot_file, patho, alpha, d_max):

    logger.info("Parsing interactome")
    interactome = data_parser.parse_interactome(interactome_file)

    logger.info("Parsing gene-to-ENSG mapping")
    ENSG2gene, gene2ENSG, Uniprot2ENSG = data_parser.parse_Uniprot(Uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, interactome, patho)

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
        Calculate GBA centrality for infertility based on the guilt-by-association paradigm.
        """
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
             patho=args.patho,
             Uniprot_file=args.Uniprot_file,
             alpha=args.alpha,
             d_max=args.d_max)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
