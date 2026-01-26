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

import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_network(network_file, weighted, directed):
    '''
    Creates an edge list representing the network

    arguments:
    - network_file: filename (with path) of network in SIF format
      (3 tab-separated columns: node1 weight/interaction_type node2)
      NOTE: weights can be floats in ]0, 1] or strings representing interaction type (eg "pp")
    - weighted: bool, if True second column must be a float in ]0,1]
    - directed: bool, if False reverse edges (with same weight) are created

    returns:
    - network: list of "edges", an edge is a tuple (source, dest, weight) where
      source and dest are ints, and weight is a float
    - node2idx: type=dict, key=node, value=unique identifier for the node, these are
      consecutive ints starting at 0
    '''
    num_nodes = 0
    num_edges = 0
    network = []
    node2idx = {}

    try:
        f = open(network_file, 'r')
    except Exception as e:
        logger.error("Opening provided SIF network file %s: %s", network_file, e)
        raise Exception("cannot open provided network file")

    # dict: key == "source&dest", value == weight for every edge that has been created
    edge2weight = {}

    for line in f:
        split_line = line.rstrip().split('\t')
        if len(split_line) != 3:
            logger.error("SIF file %s has bad line (not 3 tab-separated fields): %s",
                         network_file, line)
            raise Exception("Bad line in the network file, not 3 tab-separated fields")

        (node1, weight_or_type, node2) = split_line

        # assign unique identifiers to nodes
        if node1 not in node2idx:
            node2idx[node1] = num_nodes
            num_nodes += 1
        if node2 not in node2idx:
            node2idx[node2] = num_nodes
            num_nodes += 1

        weight = 1.0
        if weighted:
            try:
                weight = float(weight_or_type)
            except Exception:
                logger.error("SIF file %s has bad line, --weighted but weight is not a number: %s",
                             network_file, line)
                raise Exception("Bad line in the interactome file, weight is not a number")

        # create edge(s)
        edgeID = str(node2idx[node1]) + '&' + str(node2idx[node2])
        # make sure same edge wasn't seen before (with different weight)
        if edgeID in edge2weight:
            if (edge2weight[edgeID] != weight):
                logger.error("SIF file %s contains the same interaction several times with different weights, " +
                             "second occurrence is %s", network_file, line)
                raise Exception("Bad line in the interactome file")
            # else this edge was seen before with same weight, nothing to do
        else:
            # need to create this edge
            network.append((node2idx[node1], node2idx[node2], weight))
            num_edges += 1
            edge2weight[edgeID] = weight
            if not directed:
                # make sure that if reverse edge was defined, it has the same weight
                revEdgeID = str(node2idx[node2]) + '&' + str(node2idx[node1])
                if (revEdgeID in edge2weight):
                    if edge2weight[revEdgeID] != weight:
                        logger.error("SIF file %s contains the same interaction in both directions with different " +
                                     "weights, but network is undirected and weighted! Second occurrence is: %s",
                                     network_file, line)
                        raise Exception("Bad line in the interactome file")
                    # else rev edge exists and has same weight, nothing to do
                else:
                    # create reverse edge
                    network.append((node2idx[node2], node2idx[node1], weight))
                    num_edges += 1
                    edge2weight[revEdgeID] = weight
            # else network is directed, never create reverse edges

    f.close()
    logger.info("built non-redundant network with %i edges between %i nodes",
                num_edges, num_nodes)
    # sanity
    if (len(node2idx) != num_nodes):
        logger.error("sanity check failed on num_nodes!")
    if (len(network) != num_edges):
        logger.error("sanity check failed on num_edges!")
    return(network, node2idx)


def parse_seeds(seeds_file, node2idx):
    '''
    Build a vector of seeds from seeds_file

    arguments:
    - seeds_file: filename (with path) of seeds, one seed per line

    returns:
    - seeds: list of floats, one per seed, value=1 if node is a seed and 0 otherwise
    '''
    seeds = [0.0] * len(node2idx)
    num_found_seeds = 0

    try:
        f_seeds = open(seeds_file, 'r')
    except Exception as e:
        logger.error("Opening provided seeds file %s: %s",
                     seeds_file, e)
        raise Exception("cannot open provided seeds file")

    for line in f_seeds:
        seed = line.rstrip()
        if seed in node2idx:
            seeds[node2idx[seed]] = 1.0
            num_found_seeds += 1
        else:
            logger.warning("seed %s is not in the network, skipping it",
                           seed)

    f_seeds.close()

    logger.info("found %i seeds", num_found_seeds)

    return(seeds)


def scores_to_TSV(scores, node2idx):
    '''
    Print scores to stdout in TSV format, 2 columns: node score

    arguments:
    - scores: list of scores, one per node
    - node2idx: type=dict, key=node, value=index in scores
    '''

    # header
    print("NODE\tSCORE")

    for node in node2idx.keys():
        score = scores[node2idx[node]]
        print(node + "\t" + "{:.5f}".format(score))
