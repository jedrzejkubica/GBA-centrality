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


def parse_uniprot_file(uniprot_file):
    """
    Parse a TSV file from uniprot_parser.py with columns:
    - Uniprot Primary AC
    - Uniprot Secondary AC(s)
    - tax ID
    - gene name(s)

    Returns:
    - primary2secondary: dict with key=Uniprot Primary AC,
        value=list of Uniprot Secondary ACs
    - secondary2primary: dict with key=Uniprot Secondary AC, value=lsit of Uniprot Primary ACs
    """

    primary2secondary = {}
    secondary2primary = {}

    try:
        f = open(uniprot_file)
    except Exception as e:
        logger.error(f"Opening provided uniprot file {uniprot_file}: {e}")
        raise Exception(f"cannot open provided uniprot file {uniprot_file}")

    header = f.readline()

    for line in f:
        line_split = line.rstrip("\n").split("\t")

        if(len(line_split) != 4):
            logger.error(f"Uniprot file {uniprot_file} has bad line (not 4 tab-separated fields): {line}")
            raise Exception(f"Bad line in the uniprot file {uniprot_file}, not 4 tab-separated fields")
        
        primaryAC = line_split[0]
        if line_split[1] == "":
            secondaryACs = []
        else:
            secondaryACs = line_split[1].split(",")

        primary2secondary[primaryAC] = secondaryACs
        for secondaryAC in secondaryACs:
            if secondaryAC in secondary2primary:
                secondary2primary[secondaryAC].append(primaryAC)
            else:
                secondary2primary[secondaryAC] = [primaryAC]

    f.close()

    return(primary2secondary, secondary2primary)


def parse_interaction_file(interaction_file, primary2secondary, secondary2primary):
    """
    Parse a miTAB 2.5 or 2.7 file.

    For each interaction, find Uniprot ACs of interacting proteins,
    interaction detection method, pubmed and interaction type.
    Ignore lines if any of the above is missing.
    Filter interactions as follows:
    - ignore self-interactions
    - ignore some "bad" detection methods (search for "bad" below)
    - tax ID must be human (9606)
    - interaction type is used to ignore interactions or set evidence type ("1" or "2")
    sort alphabetically the two interactors (A:B and B:A are the same)

    Print to STDOUT in TSV format:
    - protein A Uniprot AC
    - protein B Uniprot AC
    - pubmed
    - evidence type
    """

    re_uniprot = re.compile(r'^uniprot(kb|/swiss-prot):([A-Z0-9-_]+)$')
    re_psimi = re.compile(r'^psi-mi:"(MI:\d+)"')  # detection method and interaction type
    re_pubmed = re.compile(r'^pubmed:(\d+)$')
    re_taxID = re.compile(r'^taxid:(\d+)')

    try:
        f = open(interaction_file, 'r')
    except Exception as e:
        logger.error(f"Opening provided interaction file {interaction_file}: {e}")
        raise Exception(f"cannot open provided interaction file {interaction_file}")

    header = f.readline()

    line_count = 0

    for line in f:
        line_count += 1
        line_split = line.rstrip().split("\t")

        # Uniprot AC of protein A should be in column 0,
        # otherwise it can be in alternatives (column 2) or in aliases (column 4);
        # for protein B, Uniprot AC should be in column 1,
        # otherwise it can be in alternatives (column 3) or in aliases (column 5)
        protein_A = ""
        potentials = [line_split[0]] + line_split[2].split("|") + line_split[4].split("|")
        for potential in potentials:
            if(re_uniprot.match(potential)):
                protein_A = re_uniprot.match(potential).group(2)  # second parenthesized group in re_uniprot
                # check if it is a valid primary AC
                if(protein_A in primary2secondary):
                    break
                if(protein_A in secondary2primary):
                    # if it is a secondary AC, keep the first primary AC associated to it
                    protein_A = secondary2primary[protein_A][0]
                    break
                else:
                    protein_A = ""
        if(protein_A == ""):
            logger.warning(f"Uniprot AC not found at line {line_count}, skipping it")
            continue

        # interactor B
        protein_B = ""
        potentials = [line_split[1]] + line_split[3].split("|") + line_split[5].split("|")
        for potential in potentials:
            if(re_uniprot.match(potential)):
                protein_B = re_uniprot.match(potential).group(2)
                # check if it is a valid primary AC
                if(protein_B in primary2secondary):
                    break
                if(protein_B in secondary2primary):
                    # if it is a secondary AC, keep the first primary AC associated to it
                    protein_B = secondary2primary[protein_B][0]
                    break
                else:
                    protein_B = ""
        if(protein_B == ""):
            logger.warning(f"Uniprot AC not found at line {line_count}, skipping it")
            continue
        
        # ignore self-interactions
        if(protein_B == protein_A):
            continue

        # interaction detection methods should be in column 6;
        # detection method cannot be "bad": MI:0254 (genetic interference) or
        # MI:0686 (unspecified method)
        method = ""
        methods_split = line_split[6].split("|")
        for met in methods_split:
            if(re_psimi.match(met)):
                method = re_psimi.match(met).group(1)
                if(method in ["MI:0254", "MI:0686"]):
                    method = ""
                    continue
                else:
                    break
        if(method == ""):
            logger.warning(f"Detection method for {protein_A}:{protein_B} not found or bad method, skipping it")
            continue

        # pubmed ID should be in column 8
        pubmed = ""
        pub_split = line_split[8].split("|")
        for pub in pub_split:
            if(re_pubmed.match(pub)):
                pubmed = re_pubmed.match(pub).group(1)
                break
        if(pubmed == ""):
            logger.warning(f"Pubmed ID not found at line {line_count}, skipping it")
            continue

        # tax ID for protein A should be in column 9, tax ID for protein B should be in column 10
        # both proteins should be human ("9606")
        taxID_A = ""
        tax_A_split = line_split[9].split("|")
        for tax in tax_A_split:
            if(re_taxID.match(tax)):
                taxID_A = re_taxID.match(tax).group(1)
                break
        if(taxID_A == ""):
            logger.warning(f"Tax ID not found at line {line_count}, skipping it")
            continue

        taxID_B = ""
        tax_B_split = line_split[10].split("|")
        for tax in tax_B_split:
            if((re_taxID.match(tax))):
                taxID_B = re_taxID.match(tax).group(1)
                break
        if(taxID_B == ""):
            logger.warning(f"Tax ID not found at line {line_count}, skipping it")
            continue
        
        # ignore non-human interactions
        if((taxID_A != "9606") or (taxID_B != "9606")):
            continue

        # interaction type should be in column 11;
        # interaction type cannot be "bad", ie MI:0403 (colocalization),
        # if interaction type is MI:0407 (direct interaction), evidence_type="1"
        # otherwise evidence_type="2"
        interaction_type = ""
        evidence_type = ""
        types_split = line_split[11].split("|")  # 19/05/2026 intact and biogrid only store one type
        for type in types_split:
            if(re_psimi.match(type)):
                interaction_type = re_psimi.match(type).group(1)
                if(interaction_type == "MI:0407"):
                    evidence_type = "1"
                    break
                elif(interaction_type != "MI:0403"):
                    evidence_type = "2"
                    break
        if(evidence_type == ""):
            logger.warning(f"Interaction type for {protein_A}:{protein_B} not found or bad type, skipping it")
            continue

        # sort alphabetically
        if(protein_B < protein_A):
            temp = protein_A
            protein_A = protein_B
            protein_B = temp

        print("\t".join([protein_A, protein_B, pubmed, evidence_type]))

    f.close()


def main(interaction_file, uniprot_file):

    logger.info("Parsing uniprot file")
    (primary2secondary, secondary2primary) = parse_uniprot_file(uniprot_file)

    logger.info("Parsing interaction file")
    parse_interaction_file(interaction_file, primary2secondary, secondary2primary)

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
        Parse a miTAB 2.5/2.7 file.
        For each interaction, find Uniprot ACs of interacting proteins,
        interaction detection method, pubmed and interaction type.
        Filter interactions based on the detection method and interaction type.
        Print to STDOUT in the TSV format:
        - protein A Uniprot AC
        - protein B Uniprot AC
        - pubmed
        - evidence type
        """)

    parser.add_argument('--interactions', required=True)
    parser.add_argument('--uniprot', required=True)

    args = parser.parse_args()

    try:
        main(interaction_file=args.interactions, uniprot_file=args.uniprot)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
