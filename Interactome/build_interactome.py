# This script was cloned from git@github.com:manojmw/grexome-TIMC-Secondary.git

import os
import sys
import argparse
import logging

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interactions(interactions_parsed_files, Uniprot2ENSG):
    """
    Parse the TSV files from interaction_parser.py with columns:
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Interaction Detection Method
    - Pubmed ID
    - Interaction Type

    Filter on Interaction Detection Method and Interaction Type:
    each interaction has at least 1 experiment,
    at least one should be proven by any binary interaction detection method.
    Removes self-loops.

    Returns a list with 5 items (in each sublist):
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Publication count
    - PubmedID(s)
    - Experiment count
    """

    PPI2PubmedID = {}
    PPI2detectionMethod = {}
    PPIs = []

    for file in interactions_parsed_files:
        f = open(file)

        for line in f:
            line_split = line.rstrip('\n').split('\t')

            detectionMethod = line_split[2]
            PubmedIDs = line_split[3]
            interactionType = line_split[4].rstrip('\n')

            # remove Interaction Detection Methods:
            # MI:0254 - genetic interference
            # MI:0686 - unspecified method
            # keep Interaction Type:
            # MI:0407 - direct interaction
            # MI:0915 - physical association
            if detectionMethod not in ['MI:0254', 'MI:0686'] and interactionType in ['MI:0407', 'MI:0915']:
                interactors = line_split[0] + '_' + line_split[1]  # proteinA + proteinB Primary Accessions

                # store PubmedIDs as a list
                # avoid duplications
                if PPI2PubmedID.get(interactors, False):
                    if PubmedIDs not in PPI2PubmedID[interactors]:
                        PPI2PubmedID[interactors].append(PubmedIDs)
                else:
                    PPI2PubmedID[interactors] = [PubmedIDs]

                # avoiding duplications
                if PPI2detectionMethod.get(interactors, False):
                    PPI2detectionMethod[interactors].append(detectionMethod)
                else:
                    PPI2detectionMethod[interactors] = [detectionMethod]

        f.close()

    for interactors in PPI2PubmedID:
        # keep PPI if at least one experiment has been proven by binary interaction method;
        # remove PPIs that have been proved only by Affintity Chromatography Technology (ACT) "MI:0004"
        if any(exp != "MI:0004" for exp in PPI2detectionMethod[interactors]):
            (proteinA, proteinB) = interactors.split('_')

            if (proteinA in Uniprot2ENSG) and (proteinB in Uniprot2ENSG):
                # remove self-loops
                if proteinA == proteinB:
                    continue
                else:
                    PubmedID = ', '.join(PPI2PubmedID[interactors])
                    PubmedID_count = len(PPI2PubmedID[interactors])
                    experiment_count = len(PPI2detectionMethod[interactors])

                    out_line = [Uniprot2ENSG[proteinA],
                                Uniprot2ENSG[proteinB],
                                str(PubmedID_count),
                                PubmedID,
                                str(experiment_count)]
                    PPIs.append(out_line)

    return PPIs


def remove_hubs(PPIs):
    """
    Parses the output from parse_interactions() list with 5 items (in each sublist):
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Publication count
    - PubmedID(s)
    - Experiment count

    Counts the degree (number of interactions) of each protein.
    Removes 1% of the highest-degree proteins.

    Returns a list with 5 items (in each sublist):
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Publication count
    - PubmedID(s)
    - Experiment count
    """

    protein_degrees = {}  # dict with key=Protein Uniprot Primary Accession, value=degree
    PPIs_noHubs = []

    for PPI in PPIs:
        proteinA = PPI[0]
        proteinB = PPI[1]

        if proteinA in protein_degrees:
            protein_degrees[proteinA] += 1
        else:
            protein_degrees[proteinA] = 1

        if proteinB in protein_degrees:
            protein_degrees[proteinB] += 1
        else:
            protein_degrees[proteinB] = 1

    # get top 1% of proteins with the highest degrees (number of interactions)
    proteins_sorted = sorted(protein_degrees.keys(), key=lambda item: protein_degrees[item], reverse=True)
    degrees_sorted = [protein_degrees[protein] for protein in proteins_sorted]

    num_elements = round(len(degrees_sorted) * 0.01)
    hubs = set(proteins_sorted[:num_elements])

    for PPI in PPIs:
        proteinA, proteinB, PubmedID_count, PubmedID, experiment_count = PPI

        if (proteinA in hubs) or (proteinB in hubs):
            continue
        else:
            out_line = [proteinA,
                        proteinB,
                        PubmedID_count,
                        PubmedID,
                        experiment_count]
            PPIs_noHubs.append(out_line)

    return PPIs_noHubs


def save_interactome(PPIs):
    """
    Parse PPIs from parse_interactions() and output from data_parser.parse_uniprot().

    Print interactome to STDOUT in SIF format:
    - ENSG of protein A
    - "pp" for "protein-protein interaction"
    - ENSG of protein B
    """
    for PPI in PPIs:
        proteinA = PPI[0]
        proteinB = PPI[1]

        out_line = (proteinA, "pp", proteinB)
        print('\t'.join(out_line))


def main(interactions_parsed_files, uniprot_file):

    logger.info("Parsing Uniprot file")
    ENSG2Gene, gene2ENSG, Uniprot2ENSG = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing interaction files")
    PPIs = parse_interactions(interactions_parsed_files, Uniprot2ENSG)

    logger.info("Removing protein hubs")
    PPIs_noHubs = remove_hubs(PPIs)

    logger.info("Printing interactome")
    save_interactome(PPIs_noHubs)

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
        Parses the output file(s) produced by interaction_parser.py
        and uniprot_parsed.tsv produced by uniprot_parser.py
        to produce a low-confidence interactome (ie, no hubs and no self-loops)
        in SIF format:
        - ENSG of Protein A
        - "pp" for "protein-protein interaction"
        - ENSG of Protein B
        """)

    parser.add_argument('--interactions_parsed_files', nargs='+', required=True)
    parser.add_argument('--uniprot_file', required=True)

    args = parser.parse_args()

    try:
        main(interactions_parsed_files=args.interactions_parsed_files,
             uniprot_file=args.uniprot_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
