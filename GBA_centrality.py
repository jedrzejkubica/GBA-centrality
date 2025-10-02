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
import ctypes

import argparse

import numpy
import networkx

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


class edge(ctypes.Structure):
    _fields_ = [('source', ctypes.c_uint),
                ('dest', ctypes.c_uint),
                ('weights', ctypes.c_float)]
    
class network(ctypes.Structure):
    _fields_ = [('nbNodes', ctypes.c_uint),
                ('nbEdges', ctypes.c_uint),
                ('edges', edge)]

class geneScores(ctypes.Structure):
    _fields_ = [('nbGenes', ctypes.c_uint),
                ('scores', ctypes.POINTER(ctypes.c_float))]


def calculate_scores(interactome, causal_genes, alpha) -> dict:
    '''
    Calculate scores for every gene in the interactome based on the proximity to known causal genes.

    arguments:
    - interactome: type=networkx.Graph
    - causal_genes: dict of causal genes with key=ENSG, value=1
    - alpha: attenuation coefficient (parameter set by user)

    returns:
    - scores: dict with key=ENSG, value=score
    '''
    so_file = "./GBA-centrality-C/gbaCentrality.so"
    gbaLibrary = ctypes.CDLL(so_file)
    # declare function signature
    gbaLibrary.gbaCentrality.argtypes = [
        ctypes.POINTER(network),
        ctypes.POINTER(geneScores),
        ctypes.c_float,
        ctypes.POINTER(geneScores)
    ]
    gbaLibrary.gbaCentrality.restype = None

    # generate adjacencyMatrix
    weights = networkx.to_numpy_array(interactome, dtype=numpy.float32)
    weights_1D = weights.flatten()
    weights_ctype = (ctypes.c_float * len(weights_1D))(*weights_1D)
    A = adjacencyMatrix(nbCols=len(interactome.nodes()),
                        weights=weights_ctype)
    
    # generate edge list
    # get node indexing instead of ENSGs
    

    # generate geneScores
    # array for genes in the interactome: 1 if causal gene, 0 otherwise
    # size=len(nodes in interactome), ordered as in interactome.nodes()
    causal_genes_vec = [0.0] * len(interactome.nodes())
    ni = 0
    for n in interactome.nodes():
        if n in causal_genes:
            causal_genes_vec[ni] = 1.0
        ni += 1

    causal_genes_ctype = (ctypes.c_float * len(causal_genes_vec))(*causal_genes_vec)
    causal = geneScores(nbGenes=len(causal_genes_vec),
                        scores=causal_genes_ctype)

    scores_vec = (ctypes.c_float * len(causal_genes_vec))()

    res = geneScores(nbGenes=len(causal_genes_vec), scores=ctypes.cast(scores_vec, ctypes.POINTER(ctypes.c_float)))

    gbaLibrary.gbaCentrality(
        ctypes.byref(A),
        ctypes.byref(causal),
        ctypes.c_float(alpha),
        ctypes.byref(res)
    )

    res_scores = [res.scores[i] for i in range(res.nbGenes)]

    scores = dict(zip(interactome.nodes(), res_scores))  # map scores to genes

    return scores


def main(interactome_file, causal_genes_file, uniprot_file, alpha, weighted, directed):

    logger.info("Parsing interactome")
    interactome = data_parser.parse_interactome(interactome_file, weighted, directed)

    sys.exit()

    logger.info("Parsing gene-to-ENSG mapping")
    ENSG2gene, gene2ENSG, uniprot2ENSG = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, interactome)

    logger.info("Calculating scores")
    scores = calculate_scores(interactome, causal_genes, alpha)

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
                        help='''filename (with path) of interactome in SIF format
                                (3 tab-separated columns: ENSG1 weight/interaction_type ENSG2), type=str
                                NOTE: second column is either weights (floats in [0, 1])
                                or one interaction type (eg "pp"); if weighted, use parameter --weighted''',
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
    parser.add_argument('--weighted',
                        help='use if graph is weighted',
                        default=False,
                        type=bool)
    parser.add_argument('--directed',
                        help='use if graph is directed',
                        default=False,
                        type=bool)

    args = parser.parse_args()

    try:
        main(interactome_file=args.interactome,
             causal_genes_file=args.causal,
             uniprot_file=args.uniprot,
             alpha=args.alpha,
             weighted=args.weighted,
             directed=args.directed)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
