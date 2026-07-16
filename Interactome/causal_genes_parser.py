############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024-2026
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
import argparse
import logging
import re

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_uniprot(uniprot_file):
    """
    Parse a TSV file from uniprot_parser.py with columns:
    - Uniprot Primary AC
    - Uniprot Secondary AC(s) (comma-separated)
    - tax ID
    - gene name

    Returns:
      - gene2uniprot: dict with key=gene name, value=Uniprot Primary AC(s)

    Note: more than one protein can be associated with a gene,
          then keeping all proteins associated with that gene as causal proteins
    """
    gene2uniprot = {}

    with open(uniprot_file, 'r') as f:
        # skip header
        header = f.readline()
        if not header.startswith("PrimaryAC\t"):
            raise Exception(f"uniprot file {uniprot_file} is headerless? expecting headers but got {header}")

        for line in f:
            line_split = line.rstrip("\n").split("\t")

            # if some records are incomplete, die
            if(len(line_split) != 4):
                raise Exception(f"Bad line in the uniprot file {uniprot_file}, not 4 tab-separated fields")

            (primaryAC, secondaryACs, taxID, gene) = line_split
            # only keep human proteins
            if taxID != "9606":
                continue

            # make sure there is a gene name
            if gene == "":
                continue
            if gene in gene2uniprot:
                gene2uniprot[gene].append(primaryAC)
            else:
                gene2uniprot[gene] = [primaryAC]

    return(gene2uniprot)


def parse_causal_genes(causal_genes_file, gene2uniprot):
    '''
    Build a list of protein Uniprot Primary AC(s) corresponding to
    causal genes from causal_genes_file

    arguments:
    - causal_genes_file: filename (with path) of known causal genes, one gene name per line
    - gene2uniprot: dict with key=gene name, value=Uniprot Primary AC(s)

    returns:
    - causal_proteins: list of Uniprot Primary AC(s)
    '''
    causal_proteins = []
    num_found_genes = 0
    num_causal_proteins = 0

    with open(causal_genes_file, 'r') as f_causal:
        for line in f_causal:
            gene = line.rstrip("\n")
            if gene in gene2uniprot:
                if len(gene2uniprot[gene]) > 1:
                    logger.warning(f"Causal gene {gene} has multiple Uniprot Primary AC(s): {gene2uniprot[gene]}")
                causal_proteins.extend(gene2uniprot[gene])
                num_causal_proteins += len(gene2uniprot[gene])
                num_found_genes += 1
            else:
                logger.warning(f"Causal gene {gene} is not a known gene in Uniprot, fix it")
                raise Exception(f"Unknown gene {gene} in Uniprot")

    logger.info(f"Found {num_causal_proteins} proteins, corresponding to {num_found_genes} causal genes")

    return(causal_proteins)


def main(uniprot_file, causal_genes_file):

    logger.info("Parsing Uniprot file")
    gene2uniprot = parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_proteins = parse_causal_genes(causal_genes_file, gene2uniprot)

    logger.info("Printing protein seeds")
    for protein in causal_proteins:
        print(protein)

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
        Parses the causal genes file and the Uniprot file to produce a list of causal proteins
        and prints causal proteins to stdout.
        """)
    parser.add_argument('--uniprot', required=True)
    parser.add_argument('--causal', required=True)

    args = parser.parse_args()

    try:
        main(args.uniprot, args.causal)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)