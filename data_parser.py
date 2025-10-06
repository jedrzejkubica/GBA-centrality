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
import re


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interactome(interactome_file, weighted, directed):
    '''
    Creates an edge list representing the interactome

    arguments:
    - interactome_file: filename (with path) of interactome in SIF format
      (3 tab-separated columns: ENSG1 weight/interaction_type ENSG2)
      NOTE: weights can be floats in ]0, 1] or strings representing interaction type (eg "pp")
    - weighted: bool, if True second column must be a float in ]0,1]
    - directed: bool, if False reverse edges (with same weight) are created

    returns:
    - interactome: list of "edges", an edge is a tuple (source, dest, weight) where
      source and dest are ints, and weight is a float
    - ENSG2idx: type=dict, key=ENSG, value=unique identifier for the ENSG, these are
      consecutive ints starting at 0
    '''
    num_nodes = 0
    num_edges = 0
    interactome = []
    ENSG2idx = {}

    try:
        f = open(interactome_file, 'r')
    except Exception as e:
        logger.error("Opening provided SIF interactome file %s: %s", interactome_file, e)
        raise Exception("cannot open provided interactome file")

    # dict: key == "source&dest", value == weight for every edge that has been created
    edge2weight = {}

    for line in f:
        split_line = line.rstrip().split('\t')
        if len(split_line) != 3:
            logger.error("SIF file %s has bad line (not 3 tab-separated fields): %s",
                         interactome_file, line)
            raise Exception("Bad line in the interactome file, not 3 tab-separated fields")

        (node1, weight_or_type, node2) = split_line

        # assign unique identifiers to nodes
        if node1 not in ENSG2idx:
            ENSG2idx[node1] = num_nodes
            num_nodes += 1
        if node2 not in ENSG2idx:
            ENSG2idx[node2] = num_nodes
            num_nodes += 1

        weight = 1.0
        if weighted:
            try:
                weight = float(weight_or_type)
            except Exception:
                logger.error("SIF file %s has bad line, --weighted but weight is not a number: %s",
                             interactome_file, line)
                raise Exception("Bad line in the interactome file, weight is not a number")

        # create edge(s)
        edgeID = str(ENSG2idx[node1]) + '&' + str(ENSG2idx[node2])
        # make sure same edge wasn't seen before (with different weight)
        if edgeID in edge2weight:
            if (edge2weight[edgeID] != weight):
                logger.error("SIF file %s contains the same interaction several times with different weights, " +
                             "second occurrence is %s", interactome_file, line)
                raise Exception("Bad line in the interactome file")
            # else this edge was seen before with same weight, nothing to do
        else:
            # need to create this edge
            interactome.append((ENSG2idx[node1], ENSG2idx[node2], weight))
            num_edges += 1
            edge2weight[edgeID] = weight
            if not directed:
                # make sure that if reverse edge was defined, it has the same weight
                revEdgeID = str(ENSG2idx[node2]) + '&' + str(ENSG2idx[node1])
                if (revEdgeID in edge2weight):
                    if edge2weight[revEdgeID] != weight:
                        logger.error("SIF file %s contains the same interaction in both directions with different " +
                                     "weights, but network is undirected and weighted! Second occurrence is: %s",
                                     interactome_file, line)
                        raise Exception("Bad line in the interactome file")
                    # else rev edge exists and has same weight, nothing to do
                else:
                    # create reverse edge
                    interactome.append((ENSG2idx[node2], ENSG2idx[node1], weight))
                    num_edges += 1
                    edge2weight[revEdgeID] = weight
            # else network is directed, never create reverse edges

    f.close()
    logger.info("built non-redundant network with %i edges between %i nodes",
                num_edges, num_nodes)
    # sanity
    if (len(ENSG2idx) != num_nodes):
        logger.error("sanity check failed on num_nodes!")
    if (len(interactome) != num_edges):
        logger.error("sanity check failed on num_edges!")
    return(interactome, ENSG2idx)


def parse_uniprot(uniprot_file):
    '''
    Parses tab-seperated Uniprot file produced by Interactome/uniprot_parser.py
    which consists of 7 columns (one record per line):
    - Uniprot Primary Accession
    - Taxonomy Identifier
    - ENST (or a comma seperated list of ENSTs)
    - ENSG (or a comma seperated list of ENSGs)
    - Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
    - GeneID (or a comma seperated list of GeneIDs)
    - Gene Name (or a comma seperated list of Gene Names)

    Returns:
      - ENSG2gene: dict with key=ENSG, value=geneName
      - gene2ENSG: dict with key=gene, value=ENSG

    Note: if more than one gene name is associated with a particular ENSG,
          then keeping the first gene name from the list
    '''
    ENSG2gene = {}
    gene2ENSG = {}

    try:
        f = open(uniprot_file, 'r')
    except Exception as e:
        logging.error("Opening provided uniprot file %s: %s", uniprot_file, e)
        raise Exception("cannot open provided Uniprot file")

    # skip header
    line = f.readline()
    if not line.startswith("Primary_AC\t"):
        logging.error("uniprot file %s is headerless? expecting headers but got %s",
                      uniprot_file, line)
        raise Exception("Uniprot file problem")

    for line in f:
        split_line = line.rstrip('\r\n').split('\t')

        # if some records are incomplete, die
        if len(split_line) != 7:
            logging.error("uniprot file %s line doesn't have 7 fields: %s",
                          uniprot_file, line)
            raise Exception("Uniprot file problem")

        (AC_primary, TaxID, ENSTs, ENSGs, AC_secondary, GeneIDs, geneNames) = split_line

        # make sure there is at least one ENSG and keep only the first one
        if ENSGs == "":
            continue
        ENSG = ENSGs.split(',')[0]

        # make sure there is at least one gene name and keep only the first one
        if geneNames == "":
            continue
        geneName = geneNames.split(',')[0]

        ENSG2gene[ENSG] = geneName
        gene2ENSG[geneName] = ENSG

    return(ENSG2gene, gene2ENSG)


def parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2idx):
    '''
    Build a vector of causal ENSGs from causal_genes_file

    arguments:
    - causal_genes_file: filename (with path) of known causal genes, one gene name per line
    - gene2ENSG: dict of all known genes, key=gene_name, value=ENSG
    - ENSG2idx: type=dict, key=ENSG, value=index in interactome

    returns:
    - causal_genes: list of floats, one per gene, value=1 if gene is causal and 0 otherwise
    '''
    causal_genes = [0.0] * len(ENSG2idx)
    num_found_genes = 0

    try:
        f_causal = open(causal_genes_file, 'r')
    except Exception as e:
        logger.error("Opening provided causal genes file %s: %s", causal_genes_file, e)
        raise Exception("cannot open provided causal genes file")

    re_causal = re.compile(r'^[a-zA-Z0-9\-_]+$')  # allow for: letters, digits, "_", "-"

    for line in f_causal:
        if re_causal.match(line):
            gene_name = line.rstrip()
            if gene_name in gene2ENSG:
                ENSG = gene2ENSG[gene_name]
                if ENSG in ENSG2idx:
                    causal_genes[ENSG2idx[ENSG]] = 1.0
                    num_found_genes += 1
                else:
                    logger.warning("causal gene %s == %s is not in interactome, skipping it",
                                   gene_name, ENSG)
            else:
                logger.warning("causal gene %s is not a known gene in gene2ENSG, skipping it",
                               gene_name)
        else:
            logger.error("Bad line in the causal genes file, doesn't look like a gene name: %s", line)
            raise Exception("Bad line in the causal genes file")

    f_causal.close()

    logger.info("found %i causal genes with known ENSG", num_found_genes)

    return(causal_genes)


def scores_to_TSV(scores, ENSG2gene, ENSG2idx):
    '''
    Print scores to stdout in TSV format, 3 columns: ENSG gene_name score

    arguments:
    - scores: list of scores, one per gene
    - ENSG2gene: dict of all known gene names, key=ENSG, value=gene_name
    - ENSG2idx: type=dict, key=ENSG, value=index in scores
    '''

    # header
    print("ENSG\tGENE\tSCORE")

    for ENSG in ENSG2idx.keys():
        gene = ENSG2gene[ENSG]
        print(ENSG + "\t" + gene + "\t" + str(scores[ENSG2idx[ENSG]]))
