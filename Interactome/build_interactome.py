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
import argparse
import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interactions(interactions_file):
    """
    Parses an interaction file with 4 columns:
    - protein A Uniprot AC
    - protein B Uniprot AC
    - pubmed
    - evidence type

    returns:
    - PPI2pubmed2method: dictionary with key: str "protein_A:protein_B",
        value: dictionary with key: str pubmedID, value: list of two ints: 
        - count of evidence type "1" (direct)
        - total count of evidences ("1" direct and "2" indirect)
    """
    PPI2pubmed2method = {}
    
    with open(interactions_file, 'r') as f:
        for line in f:
            line_split = line.rstrip().split('\t')

            if(len(line_split) != 4):
                logger.error(f"Interactions file {interactions_file} has bad line (not 4 tab-separated fields): {line}")
                raise Exception(f"Bad line in the interactions file {interactions_file}, not 4 tab-separated fields")

            PPI = line_split[0] + ':' + line_split[1]  # str "protein_A:protein_B"
            pubmed = line_split[2]
            evidence_type = line_split[3]

            if PPI not in PPI2pubmed2method:
                PPI2pubmed2method[PPI] = {}
            if pubmed not in PPI2pubmed2method[PPI]:
                PPI2pubmed2method[PPI][pubmed] = [0, 0]
            if evidence_type == "1":
                PPI2pubmed2method[PPI][pubmed][0] += 1
                PPI2pubmed2method[PPI][pubmed][1] += 1
            elif evidence_type == "2":
                PPI2pubmed2method[PPI][pubmed][1] += 1

    return(PPI2pubmed2method)


def main(interactions_parsed_files, n_evidence, n_direct):
    # merge interactions from all interaction files (eg. BioGRID, IntAct) into one dictionary:
    # key: str "protein_A:protein_B"
    # value: dictionary with key: str pubmedID, value: list of two ints:
    # - max count of evidence type "1" (direct) from any file
    # - max total count of evidences ("1" direct and "2" indirect) from any file
    PPI2pubmed2method_merged = {}

    # we have multiple files, one line is an interaction with publication
    # and evidence type ("1" or "2"),
    # interactions can be redundant between files,
    # one publication (ie. pubmedID) can report more than one interaction
    for file in interactions_parsed_files:
        logger.info(f"Parsing {file}")
        PPI2pubmed2method = parse_interactions(file)

        for PPI in PPI2pubmed2method:
            if PPI not in PPI2pubmed2method_merged:
                PPI2pubmed2method_merged[PPI] = PPI2pubmed2method[PPI].copy()
            else:
                for pubmed in PPI2pubmed2method[PPI]:
                    if pubmed not in PPI2pubmed2method_merged[PPI]:
                        PPI2pubmed2method_merged[PPI][pubmed] = PPI2pubmed2method[PPI][pubmed].copy()
                    else:
                        # keep max of direct and max of total
                        if PPI2pubmed2method[PPI][pubmed][0] > PPI2pubmed2method_merged[PPI][pubmed][0]:
                            PPI2pubmed2method_merged[PPI][pubmed][0] = PPI2pubmed2method[PPI][pubmed][0]
                        if PPI2pubmed2method[PPI][pubmed][1] > PPI2pubmed2method_merged[PPI][pubmed][1]:
                            PPI2pubmed2method_merged[PPI][pubmed][1] = PPI2pubmed2method[PPI][pubmed][1]
    
    logger.info(f"Filtering on evidence (n_evidence>={n_evidence}, n_direct>={n_direct})")
    for PPI in PPI2pubmed2method_merged:
        # sum evidences for each interaction, all publications (ie. pubmedIDs) combined
        evidence_sum = [0, 0]  # sum of direct and sum of total
        for pubmed in PPI2pubmed2method_merged[PPI]:
            evidence_sum[0] += PPI2pubmed2method_merged[PPI][pubmed][0]
            evidence_sum[1] += PPI2pubmed2method_merged[PPI][pubmed][1]
        
        # keep interactions with at least N total evidences including N direct
        if (evidence_sum[1] >= n_evidence) and (evidence_sum[0] >= n_direct):
            (protein_A, protein_B) = PPI.split(':')
            print('\t'.join([protein_A, "pp", protein_B]))

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
        Builds and prints an interactome to STDOUT in a SIF-like format:
        - protein A
        - "pp" (for "protein-protein interaction")
        - protein B
        """)

    parser.add_argument('--interactions', nargs='+', required=True)
    parser.add_argument('--n_evidence', required=False, default=2, type=int)
    parser.add_argument('--n_direct', required=False, default=1, type=int)

    args = parser.parse_args()

    try:
        main(interactions_parsed_files=args.interactions,
             n_evidence=args.n_evidence,
             n_direct=args.n_direct)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)