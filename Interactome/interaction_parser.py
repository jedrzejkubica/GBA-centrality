# This script is based on git@github.com:manojmw/grexome-TIMC-Secondary.git

############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2026
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
import re
import argparse
import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interaction_file(interaction_file):
    """
    Parse a miTAB 2.5 or 2.7 file.

    For each interaction, find Uniprot IDs of interacting proteins,
    interaction detection method, pubmed and interaction type.
    Ignore lines if any protein and/or pubmed is missing.
    Filter interactions as follows:
    - ignore self-interactions
    - tax ID must be human (9606)
    - detection method cannot be MI:0254 (genetic interference), MI:0686 (unspecified method)
        or MI:0004 (affintity chromatography)
    - interaction type must be MI:0407 (direct interaction) or MI:0915 (physical association)

    Print to STDOUT in TSV format:
    - protein A Uniprot ID
    - protein B Uniprot ID
    - interaction detection method
    - pubmed
    - interaction type
    """

    interactions = []

    re_uniprot = re.compile(r'^uniprot(kb|/swiss-prot):([A-Z0-9-_]+)$')
    re_psimi = re.compile(r'^psi-mi:"(MI:\d+)"')  # detection method and interaction type
    re_pubmed = re.compile(r'^pubmed:(\d+)$')
    re_taxID = re.compile(r'^taxid:(\d+)')

    try:
        f = open(interaction_file, 'r')
    except Exception as e:
        logger.error("Opening provided interaction file %s: %s", interaction_file, e)
        raise Exception("cannot open provided interaction file")

    header = f.readline()

    line_count = 0

    for line in f:
        line_count += 1
        line_split = line.rstrip().split("\t")

        # Uniprot ID of protein A should be in column 0,
        # otherwise it can be in alternatives (column 2) or in aliases (column 4);
        # for protein B, Uniprot ID should be in column 1,
        # otherwise it can be in alternatives (column 3) or in aliases (column 5)
        protein_A = ""
        if re_uniprot.match(line_split[0]):
            protein_A = re_uniprot.match(line_split[0]).group(2)  # second parenthesized group in re_uniprot
        else:
            alt_split = line_split[2].split("|")
            for alt in alt_split:
                if re_uniprot.match(alt):
                    protein_A = re_uniprot.match(alt).group(2)
                    break
            if protein_A == "":
                alias_split = line_split[4].split("|")
                for alias in alias_split:
                    if re_uniprot.match(alias):
                        protein_A = re_uniprot.match(alias).group(2)
                        break
        if protein_A == "":
            logger.warning(f"Uniprot ID not found at line {line_count}, skipping it")
            continue

        # interactor B
        protein_B = ""
        if re_uniprot.match(line_split[1]):
            protein_B = re_uniprot.match(line_split[1]).group(2)
        else:
            alt_split = line_split[3].split("|")
            for alt in alt_split:
                if re_uniprot.match(alt):
                    protein_B = re_uniprot.match(alt).group(2)
                    break
            if protein_B == "":
                alias_split = line_split[5].split("|")
                for alias in alias_split:
                    if re_uniprot.match(alias):
                        protein_B = re_uniprot.match(alias).group(2)
                        break
        if protein_B == "":
            logger.warning(f"Uniprot ID not found at line {line_count}, skipping it")
            continue
        
        # ignore self-interactions
        if protein_B == protein_A:
            continue

        # interaction detection methods should be in column 6
        method = ""
        methods_split = line_split[6].split("|")
        for met in methods_split:
            if re_psimi.match(met):
                method = re_psimi.match(met).group(1)
                break
        if method == "":
            logger.warning(f"Detection method for {protein_A}:{protein_B} not found, skipping it")
            continue
        if method in ["MI:0254", "MI:0686", "MI:0004"]:
            logger.warning(f"{protein_A}:{protein_B} detected by {method}, skipping it")
            continue

        # pubmed ID should be in column 8
        pubmed = ""
        pub_split = line_split[8].split("|")
        for pub in pub_split:
            if re_pubmed.match(pub):
                pubmed = re_pubmed.match(pub).group(1)
                break
        if pubmed == "":
            logger.warning(f"Pubmed ID not found at line {line_count}, skipping it")
            continue

        # tax ID for protein A should be in column 9, tax ID for protein B should be in column 10
        taxID_A = ""
        tax_A_split = line_split[9].split("|")
        for tax in tax_A_split:
            if (re_taxID.match(tax)):
                taxID_A = re_taxID.match(tax).group(1)
                break
        if taxID_A == "":
            logger.warning(f"Tax ID not found at line {line_count}, skipping it")
            continue

        taxID_B = ""
        tax_B_split = line_split[10].split("|")
        for tax in tax_B_split:
            if (re_taxID.match(tax)):
                taxID_B = re_taxID.match(tax).group(1)
        if taxID_B == "":
            logger.warning(f"Tax ID not found at line {line_count}, skipping it")
            continue
        
        # ignore non-human interactions
        if (taxID_A != "9606") or (taxID_B != "9606"):
            continue

        # interaction type should be in column 11
        interaction_type = ""
        types_split = line_split[11].split("|")
        for type in types_split:
            if re_psimi.match(type):
                interaction_type = re_psimi.match(type).group(1)
                break
        if interaction_type == "":
            logger.warning(f"Interaction type for {protein_A}:{protein_B} not found, skipping it")
            continue
        if interaction_type not in ["MI:0407", "MI:0915"]:
            # logger.warning(f"{protein_A}:{protein_B} is a {interaction_type}, skipping it")
            continue

        interactions.append((protein_A, protein_B, method, pubmed, interaction_type))
    f.close()

    for interaction in interactions:
        interaction_out_line = "\t".join(interaction) + "\n"
        print(interaction_out_line)


def main(interaction_file):

    logger.info("Parsing interaction file")
    parse_interaction_file(interaction_file)

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
        description="""
        Parse a miTAB 2.5/2.7 file, and print to STDOUT in the TSV format:
        - protein A Uniprot ID
        - protein B Uniprot ID
        - interaction detection method
        - pubmed
        - interaction type
        """)

    parser.add_argument('--interactions', required=True)

    args = parser.parse_args()

    try:
        main(interaction_file=args.interactions)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
