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


def parse_hgnc(hgnc_file):
    """
    Parse a TSV file from HGNC.
    For each entry, find a gene-to-protein mapping;
    if multiple Primary ACs are found for a gene symbol, keep all of them in a list

    Returns:
      - gene2uniprot: dict with key=gene symbol,
      value=list of Primary ACs
    """
    gene2uniprot = {}

    re_gene = re.compile(r'^([a-zA-Z0-9\-_]+)$')  # allow for: letters, digits, "_", "-"
    re_uniprot = re.compile(r'^([A-Z0-9-_]+)$')

    with open(hgnc_file, 'r') as f:
        # skip header
        header = f.readline()
        if not header.startswith("hgnc_id\t"):
            raise Exception(f"HGNC file {hgnc_file} is headerless? expecting headers but got {header}")

        for line in f:
            gene_symbols = []
            primary_ACs = []
    
            line_split = line.rstrip("\n").split("\t")

            # gene symbols should be in column 1,
            # optionally can have aliases (column 8) or previous symbols (column 10)
            if line_split[1] == "":
                logger.error(f"HGNC file {hgnc_file} has no gene symbol in line: {line}")
                raise Exception(f"Bad line in the HGNC file {hgnc_file}, no gene symbol")
            gene_symbols.append(line_split[1].rstrip())
            if line_split[8] != "":
                aliases = line_split[8].strip('"').split("|")
                for alias in aliases:
                    if re_gene.match(alias):
                        gene_symbols.append(re_gene.match(alias).group(1))
                    else:
                        logger.warning(f"Bad gene symbol {alias}, skipping it")

            if line_split[10] != "":
                previous_symbols = line_split[10].strip('"').split("|")
                for symbol in previous_symbols:
                    if re_gene.match(symbol):
                        gene_symbols.append(re_gene.match(symbol).group(1))
                    else:
                        logger.warning(f"Bad gene symbol {symbol}, skipping it")

            # Primary ACs should be in column 25
            if line_split[25] == "":
                logger.warning(f"No Primary AC found for gene {gene_symbols[0]}, skipping it")
                continue
            potential_primary_ACs = line_split[25].strip('"').split("|")
            for potential_AC in potential_primary_ACs:
                if re_uniprot.match(potential_AC):
                    primary_ACs.append(re_uniprot.match(potential_AC).group(1))
                else:
                    logger.warning(f"Bad Primary AC {potential_AC} for gene {gene_symbols[0]}, skipping it")

            for gene in gene_symbols:
                if gene in gene2uniprot:
                    gene2uniprot[gene].extend(primary_ACs)
                else:
                    gene2uniprot[gene] = primary_ACs

    return(gene2uniprot)


def parse_causal_genes(causal_genes_file, gene2uniprot):
    '''
    Build a list of protein Uniprot Primary ACs corresponding to
    causal genes from causal_genes_file

    arguments:
    - causal_genes_file: filename (with path) of known causal genes, one gene name per line
    - gene2uniprot: dict with key=gene symbol,
      value=Uniprot AC or a list of Uniprot ACs

    returns:
    - causal_proteins: list of Uniprot Primary ACs
    '''
    causal_proteins = []
    num_found_genes = 0
    num_causal_proteins = 0
    
    re_gene = re.compile(r'^[a-zA-Z0-9\-_]+$')  # allow for: letters, digits, "_", "-"

    with open(causal_genes_file, 'r') as f_causal:
        for line in f_causal:
            gene = line.rstrip("\n")
            if re_gene.match(gene):
                if gene in gene2uniprot:
                    causal_proteins.extend(gene2uniprot[gene])
                    num_causal_proteins += len(gene2uniprot[gene])
                    if len(gene2uniprot[gene]) > 1:
                        logger.warning(f"Multiple ACs found for gene {gene}: {gene2uniprot[gene]}")
                    num_found_genes += 1
                else:
                    logger.warning(f"Causal gene {gene} is not a known gene in HGNC, skipping it")
            else:
                logger.error(f"Bad line in causal genes file {causal_genes_file}, doesn't look like a gene name: {line}")
                raise Exception(f"Bad line in the causal genes file {causal_genes_file}")

    logger.info(f"Found {num_causal_proteins} proteins, corresponding to {num_found_genes} causal genes")

    return(causal_proteins)


def main(hgnc_file, causal_genes_file):

    logger.info("Parsing HGNC file")
    gene2uniprot = parse_hgnc(hgnc_file)

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
        Parses the causal genes file and the HGNC file to produce a list of causal proteins
        and prints causal proteins to stdout.
        """)
    parser.add_argument('--hgnc', required=True)
    parser.add_argument('--causal', required=True)

    args = parser.parse_args()

    try:
        main(args.hgnc, args.causal)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)