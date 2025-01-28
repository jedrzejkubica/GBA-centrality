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


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interactome(interactome_file) -> networkx.Graph:
    '''
    Creates a networkx.Graph representing the interactome

    arguments:
    - interactome_file: filename (with path) of interactome in SIF format (ie
      3 tab-separated columns: ENSG1 pp ENSG2), type=str

    returns:
    - interactome: type=networkx.Graph
    '''
    interactome = networkx.Graph()
    num_edges = 0

    try:
        f = open(interactome_file, 'r')
    except Exception as e:
        logger.error("Opening provided SIF interactome file %s: %s", interactome_file, e)
        raise Exception("cannot open provided interactome file")

    for line in f:
        split_line = line.rstrip().split('\t')
        if len(split_line) != 3:
            logger.error("SIF file %s has bad line (not 3 tab-separated fields): %s",
                         interactome_file, line)
            raise Exception("Bad line in the interactome file")

        gene1, pp, gene2 = split_line

        # exclude self-interactions
        if gene1 == gene2:
            continue
        # else: ppopulate structures
        interactome.add_edge(gene1, gene2)
        num_edges += 1

    logger.info("built non-redundant network from %i non-self interactions, " +
                "resulting in %i edges between %i nodes",
                num_edges, len(interactome.edges()), len(interactome.nodes))
    return interactome


def parse_Uniprot(Uniprot_file):
    '''
    Parses tab-seperated Uniprot file produced by Uniprot_parser.py
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
      - Uniprot2ENSG: dict with key=Primary accession, value=ENSG

    Note: if more than one gene name is associated with a particular ENSG,
          then keeping the first gene name from the list
    '''
    ENSG2gene = {}
    gene2ENSG = {}
    Uniprot2ENSG = {}

    try:
        f = open(Uniprot_file, 'r')
    except Exception as e:
        logging.error("Opening provided gene2ENSG file %s: %s", Uniprot_file, e)
        raise Exception("cannot open provided Uniprot file")

    # skip header
    next(f)

    for line in f:
        split_line = line.rstrip().split('\t')

        # if some records are incomplete, then skip them
        if len(split_line) != 7:
            continue

        AC_primary, TaxID, ENSTs, ENSGs, AC_secondary, GeneIDs, geneNames = split_line

        # keep only the first ENSG (if exists in the file)
        if ENSGs == "":
            continue
        ENSGs_list = ENSGs.rstrip().split(',')
        ENSG = ENSGs_list[0]

        # if more than one gene name is associated with a particular ENSG,
        # then keep only the first gene name from the list
        geneNames_list = geneNames.rstrip().split(',')
        geneName = geneNames_list[0]

        ENSG2gene[ENSG] = geneName
        gene2ENSG[geneName] = ENSG
        Uniprot2ENSG[AC_primary] = ENSG

    return ENSG2gene, gene2ENSG, Uniprot2ENSG


def parse_causal_genes(causal_genes_file, gene2ENSG, interactome, patho) -> dict:
    '''
    Build a dict of causal ENSGs for patho

    arguments:
    - causal_genes_file: filename (with path) of known causal genes TSV file
      with 2 columns: gene_name pathologyID, type=str
    - gene2ENSG: dict of all known genes, key=gene_name, value=ENSG
    - interactome: networkx.Graph
    - pathologyID of interest, causal genes for other pathologyIDs are ignored

    returns:
    - causal_genes: dict of all causal genes present in interactome, with key=ENSG, value=1
    '''
    causal_genes = {}

    try:
        f_causal = open(causal_genes_file, 'r')
    except Exception as e:
        logger.error("Opening provided causal genes file %s: %s", causal_genes_file, e)
        raise Exception("cannot open provided causal genes file")

    for line in f_causal:
        split_line = line.rstrip().split('\t')
        if len(split_line) != 2:
            logger.error("causal_genes file %s has bad line (not 2 tab-separated fields): %s",
                         causal_genes_file, line)
            raise Exception("Bad line in the causal_genes file")
        gene_name, pathology = split_line

        if pathology != patho:
            continue
        if gene_name in gene2ENSG:
            ENSG = gene2ENSG[gene_name]
            if ENSG in interactome:
                causal_genes[ENSG] = 1
            else:
                logger.warning("causal gene %s == %s is not in interactome, skipping it",
                               gene_name, ENSG)

    f_causal.close()

    logger.info("found %i causal genes with known ENSG for pathology %s",
                len(causal_genes), patho)

    return causal_genes


def scores_to_TSV(scores, ENSG2gene):
    '''
    Print scores to stdout in TSV format, 3 columns: ENSG gene_name score

    arguments:
    - scores: dict with key=ENSG, value=score
    - ENSG2gene: dict of all known gene names, key=ENSG, value=gene_name
    '''

    # header
    print("ENSG\tGENE\tSCORE")

    for (ENSG, score) in sorted(scores.items()):
        # GENE defaults to "" if we don't know the gene name of ENSG
        gene = ""
        if ENSG in ENSG2gene:
            gene = ENSG2gene[ENSG]
        print(ENSG + "\t" + gene + "\t" + str(score))
