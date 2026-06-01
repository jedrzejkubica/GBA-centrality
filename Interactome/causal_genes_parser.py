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
    - gene name(s) (comma-separated)

    Returns:
      - gene2uniprot: dict with key=gene name, value=Uniprot Primary AC

    Note: if more than one gene name is associated with a particular protein,
          then keeping the first gene name from the list
    """
    gene2uniprot = {}

    with open(uniprot_file, 'r') as f:
        # skip header
        header = f.readline()
        if not header.startswith("PrimaryAC\t"):
            raise Exception("uniprot file %s is headerless? expecting headers but got %s",
                        uniprot_file, header)

        for line in f:
            line_split = line.rstrip("\n").split("\t")

            # if some records are incomplete, die
            if(len(line_split) != 4):
                raise Exception("Bad line in the uniprot file, not 4 tab-separated fields")

            (primaryAC, secondaryACs, taxID, gene_names) = line_split

            # make sure there is at least one gene name and keep only the first one
            if gene_names == "":
                continue
            gene = gene_names.split(',')[0]
            gene2uniprot[gene] = primaryAC

    return(gene2uniprot)


def parse_causal_genes(causal_genes_file, gene2uniprot):
    '''
    Build a list of protein Uniprot Primary AC corresponding to
    causal gene names from causal_genes_file

    arguments:
    - causal_genes_file: filename (with path) of known causal genes, one gene name per line
    - gene2uniprot: dict mapping gene name to Uniprot Primary AC

    returns:
    - causal_proteins: list of Uniprot Primary accession for causal genes
    '''
    causal_proteins = []
    num_found_genes = 0

    with open(causal_genes_file, 'r') as f_causal:
        re_causal = re.compile(r'^[a-zA-Z0-9\-_]+$')  # allow for: letters, digits, "_", "-"

        for line in f_causal:
            gene_name = line.rstrip()
            if re_causal.match(gene_name):
                if gene_name in gene2uniprot:
                    uniprot_ac = gene2uniprot[gene_name]
                    causal_proteins.append(uniprot_ac)
                    num_found_genes += 1
                else:
                    logger.warning(f"causal gene {gene_name} is not a known gene in gene2uniprot, skipping it")
            else:
                logger.error(f"Bad line in the causal genes file, doesn't look like a gene name: {line}")
                raise Exception("Bad line in the causal genes file")

    logger.info(f"found {num_found_genes} causal genes")

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
        Parses the causal genes file and uniprot_parsed.tsv produced by uniprot_parser.py,
        print causal proteins to stdout
        """)
    parser.add_argument('--uniprot', required=True)
    parser.add_argument('--causal', required=True)

    args = parser.parse_args()

    try:
        main(args.uniprot, args.causal)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)