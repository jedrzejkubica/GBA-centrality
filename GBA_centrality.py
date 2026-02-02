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
SCORETYPE = ctypes.c_float


# following classes must match those in GBA-centrality-C/network.h and
# GBA-centrality-C/scores.h
class Edge(ctypes.Structure):
    _fields_ = [('source', ctypes.c_uint),
                ('dest', ctypes.c_uint),
                ('weight', ctypes.c_float)]


class Network(ctypes.Structure):
    _fields_ = [('nbNodes', ctypes.c_ulong),
                ('nbEdges', ctypes.c_ulong),
                ('edges', ctypes.POINTER(Edge))]


class nodeScores(ctypes.Structure):
    _fields_ = [('nbSeeds', ctypes.c_size_t),
                ('scores', ctypes.POINTER(SCORETYPE))]


def calculate_scores(network, node2idx, seeds, alpha, pathToCode, threads):
    '''
    Calculate scores for every node in the network based on the proximity to the seeds.

    arguments:
    - nodes: list of "edges", an edge is a tuple (source, dest, weight) where
      source and dest are ints, and weight is a float
    - node2idx: type=dict, key=node, value=unique identifier for the node, these are
      consecutive ints starting at 0
    - seeds: list of floats of length num_nodes, value=1 if node in seeds and 0 otherwise
    - alpha: attenuation coefficient (parameter set by user)
    - threads: number of threads to use, 0 to use all available cores

    returns:
    - scores: list of floats of length num_nodes, value=score in the same order as seeds
    '''
    if threads:
        os.environ['OMP_NUM_THREADS'] = str(threads)
    so_file = pathToCode + "/GBA-centrality-C/gbaCentrality.so"
    gbaLibrary = ctypes.CDLL(so_file)
    # declare function signature
    gbaLibrary.gbaCentrality.argtypes = [
        ctypes.POINTER(Network),
        ctypes.POINTER(nodeScores),
        ctypes.c_float,
        ctypes.POINTER(nodeScores)
    ]
    gbaLibrary.gbaCentrality.restype = None

    # generate ctypes edges
    edgesType = Edge * len(network)
    edges = edgesType()
    edgeIndex = 0
    for interaction in network:
        (source, dest, weight) = interaction
        e = Edge(ctypes.c_uint(source),
                 ctypes.c_uint(dest),
                 ctypes.c_float(weight))
        edges[edgeIndex] = e
        edgeIndex += 1

    # generate ctypes network
    N = Network(ctypes.c_ulong(len(node2idx)),
                ctypes.c_ulong(len(network)),
                edges)

    # generate ctypes seeds and scores
    seedsType = SCORETYPE * len(node2idx)
    seedsData = seedsType()
    scoresData = seedsType()

    # fill seedsData
    seedIndex = 0
    for nodeIsSeed in seeds:
        seedsData[seedIndex] = SCORETYPE(nodeIsSeed)
        seedIndex += 1

    seeds_vector = nodeScores(ctypes.c_size_t(len(node2idx)),
                        seedsData)
    scores = nodeScores(ctypes.c_size_t(len(node2idx)),
                        scoresData)

    gbaLibrary.gbaCentrality(
        ctypes.byref(N),
        ctypes.byref(seeds_vector),
        ctypes.c_float(alpha),
        ctypes.byref(scores)
    )

    scoresList = [scores.scores[i] for i in range(scores.nbSeeds)]
    return(scoresList)


def main(network_file, seeds_file, alpha, weighted, directed, pathToCode, threads):

    logger.info("Parsing network")
    (network, node2idx) = data_parser.parse_network(network_file, weighted, directed)

    logger.info("Parsing seeds")
    (seeds, seeds_vector) = data_parser.parse_seeds(seeds_file, node2idx)
    if len(seeds) > 0:
        logger.info("Calculating scores")
        scores = calculate_scores(network, node2idx, seeds_vector, alpha, pathToCode, threads)

        logger.info("Printing scores")
        data_parser.scores_to_TSV(scores, node2idx)
    else:
        logger.info("No seeds, nothing to do")

    logger.info("Done!")


if __name__ == "__main__":
    (pathToCode, script_name) = os.path.split(os.path.realpath(sys.argv[0]))
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="""
        GBA centrality is a new network propagation algorithm
        based on non-backtracking walks and in-degree normalization.
        The method assigns scores to nodes in the network that represent
        their likelihood of being in proximity to a given list of nodes of interest.
        """
    )

    parser.add_argument('--network',
                        help='''filename (with path) of network in a SIF-like format
                                (3 tab-separated columns: node1 weight/interaction_type node2), type=str
                                NOTE: second column is either weights (floats in ]0, 1])
                                or a single interaction type (eg "pp"); if weighted, use parameter --weighted''',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--seeds',
                        help='TXT file (without a header) with 1 seed per line',
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
    parser.add_argument('--threads',
                        help='number of parallel threads to run, default=0 to use all available cores',
                        default=0,
                        type=int)
 
    args = parser.parse_args()

    try:
        main(args.network, args.seeds, args.alpha, args.weighted,
             args.directed, pathToCode, args.threads)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
