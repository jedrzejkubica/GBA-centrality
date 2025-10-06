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

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

# SCORETYPE has to be identical to SCORETYPE in GBA-centrality-C/scores.h
SCORETYPE = ctypes.c_double


# following classes must match those in GBA-centrality-C/network.h and
# GBA-centrality-C/scores.h
class Edge(ctypes.Structure):
    _fields_ = [('source', ctypes.c_uint),
                ('dest', ctypes.c_uint),
                ('weight', ctypes.c_float)]


class Network(ctypes.Structure):
    _fields_ = [('nbNodes', ctypes.c_uint),
                ('nbEdges', ctypes.c_uint),
                ('edges', ctypes.POINTER(Edge))]


class GeneScores(ctypes.Structure):
    _fields_ = [('nbGenes', ctypes.c_uint),
                ('scores', ctypes.POINTER(SCORETYPE))]


def calculate_scores(interactome, ENSG2idx, causal_genes, alpha):
    '''
    Calculate scores for every gene in the interactome based on the proximity to known causal genes.

    arguments:
    - interactome: list of "edges", an edge is a tuple (source, dest, weight) where
      source and dest are ints, and weight is a float
    - ENSG2idx: type=dict, key=ENSG, value=unique identifier for the ENSG, these are
      consecutive ints starting at 0
    - causal_genes: list of floats of length num_genes, value=1 if gene is causal and 0 otherwise
    - alpha: attenuation coefficient (parameter set by user)

    returns:
    - scores: list of floats of length num_genes, value=score in the same order as causal_genes
    '''
    so_file = "./GBA-centrality-C/gbaCentrality.so"
    gbaLibrary = ctypes.CDLL(so_file)
    # declare function signature
    gbaLibrary.gbaCentrality.argtypes = [
        ctypes.POINTER(Network),
        ctypes.POINTER(GeneScores),
        ctypes.c_float,
        ctypes.POINTER(GeneScores)
    ]
    gbaLibrary.gbaCentrality.restype = None

    # generate ctypes edges
    edgesType = Edge * len(interactome)
    edges = edgesType()
    edgeIndex = 0
    for interaction in interactome:
        (source, dest, weight) = interaction
        e = Edge(ctypes.c_uint(source),
                 ctypes.c_uint(dest),
                 ctypes.c_float(weight))
        edges[edgeIndex] = e
        edgeIndex += 1

    # generate ctypes network
    N = Network(ctypes.c_uint(len(ENSG2idx)),
                ctypes.c_uint(len(interactome)),
                ctypes.pointer(edges))

    # generate ctypes causal genes and scores
    causalGenesType = SCORETYPE * len(ENSG2idx)
    causalData = causalGenesType()
    scoresData = causalGenesType()

    # fill causalData
    causalIndex = 0
    for geneIsCausal in causal_genes:
        causalData[causalIndex] = SCORETYPE(geneIsCausal)
        causalIndex += 1

    causal = GeneScores(ctypes.c_uint(len(ENSG2idx)),
                        ctypes.pointer(causalData))
    scores = GeneScores(ctypes.c_uint(len(ENSG2idx)),
                        ctypes.pointer(scoresData))

    gbaLibrary.gbaCentrality(
        ctypes.byref(N),
        ctypes.byref(causal),
        ctypes.c_float(alpha),
        ctypes.byref(scores)
    )

    scoresList = [scores.scores[i] for i in range(scores.nbGenes)]
    return(scoresList)


def main(interactome_file, causal_genes_file, uniprot_file, alpha, weighted, directed):

    logger.info("Parsing interactome")
    (interactome, ENSG2idx) = data_parser.parse_interactome(interactome_file, weighted, directed)

    logger.info("Parsing gene-to-ENSG mapping")
    (ENSG2gene, gene2ENSG) = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2idx)

    logger.info("Calculating scores")
    scores = calculate_scores(interactome, ENSG2idx, causal_genes, alpha)

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
                        action='store_true')  # if present, set the value to True; otherwise False
    parser.add_argument('--directed',
                        help='use if graph is directed',
                        action='store_true')

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
