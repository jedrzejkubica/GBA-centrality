# This script was cloned from git@github.com:manojmw/grexome-TIMC-Secondary.git

import os
import sys
import argparse
import logging

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import data_parser


def parse_interactions(interactions_parsed_files):
    """
    Parse the TSV files from interaction_parser.py with columns:
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Interaction Detection Method
    - Pubmed ID
    - Interaction Type

    Filter on Interaction Detection Method and Interaction Type:
    each interaction has at least 2 experiments,
    at least one should be proven by any binary interaction detection method

    Returns a list with 5 items (in each sublist):
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Publication count
    - PubmedID(s)
    - Experiment count
    """

    logging.info("Starting to run")

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

            # keep Interaction Detection Method:
            # MI:0096 - pull down
            # MI:0254 - genetic interference
            # MI:0686 - unspecified method
            # keep Interaction Type:
            # MI:0407 - direct interaction
            # MI:0915 - physical association
            if detectionMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and interactionType in ['MI:0407', 'MI:0915']:
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

            PubmedID = ', '.join(PPI2PubmedID[interactors])
            PubmedID_count = len(PPI2PubmedID[interactors])
            experiment_count = len(PPI2detectionMethod[interactors])

            if experiment_count >= 2:
                out_line = [proteinA, proteinB, str(PubmedID_count), PubmedID, str(experiment_count)]
                PPIs.append(out_line)

    return PPIs


def get_PPI_interactors(PPIs):
    """
    Parses the PPI list returned by parse_interactions():
    Each sublist must contain:
    - Uniprot Primary Accession of protein A
    - Uniprot Primary Accession of protein B

    Returns:
    - proteinA2interactors: dict with key=proteinA Uniprot Accession, value=list of proteinA interactors
    - proteinB2interactors: dict with key=proteinB Uniprot Accession, value=list of proteinB interactors
    """

    proteinA2interactors = {}
    proteinB2interactors = {}

    for PPI in PPIs:
        # skip self-interaction
        if PPI[0] == PPI[1]:
            continue
        else:
            if proteinA2interactors.get(PPI[0], False):
                proteinA2interactors[PPI[0]].append(PPI[1])
            else:
                proteinA2interactors[PPI[0]] = [PPI[1]]

            if proteinB2interactors.get(PPI[1], False):
                proteinB2interactors[PPI[1]].append(PPI[0])
            else:
                proteinB2interactors[PPI[1]] = [PPI[0]]

    return proteinA2interactors, proteinB2interactors


def get_protein_hubs(proteinA2interactors, proteinB2interactors):
    """
    Parses proteinA2interactors and proteinB2interactors from get_PPI_interactors().

    Checks if a protein is a hub (if it has >120 interactions).

    Returns a dictionary:
    - protein_hubs: dict with key=Uniprot Accession, value=1 if hub
    """

    proteinHubs = {}

    for protein in proteinA2interactors:
        intCountProteinA = len(proteinA2interactors[protein])
        if protein in proteinB2interactors:
            intCountProteinB = 0
            # avoid duplicated interaction count for protein A and protein B
            for interactor in proteinB2interactors[protein]:
                if interactor not in proteinA2interactors[protein]:
                    intCountProteinB += 1
            intCount = intCountProteinA + intCountProteinB
        else:
            intCount = intCountProteinA

        if intCount > 120:
            proteinHubs[protein] = 1

    for protein in proteinB2interactors:
        if protein not in proteinA2interactors:
            if len(proteinB2interactors[protein]) > 120:
                proteinHubs[protein] = 1

    return proteinHubs


def Interactome_Uniprot2ENSG(args):
    """
    Parse PPIs from parse_interactions() and output from data_parser.parse_uniprot().

    Print interactome to STDOUT in SIF format:
    - ENSG of protein A
    - "pp" for "protein-protein interaction"
    - ENSG of protein B
    """

    PPIs = parse_interactions(args.interactions_parsed_files)
    ENSG2Gene, gene2ENSG, Uniprot2ENSG = data_parser.parse_uniprot(args.uniprot_file)

    (proteinA2interactors, proteinB2interactors) = get_PPI_interactors(PPIs)
    proteinHubs = get_protein_hubs(proteinA2interactors, proteinB2interactors)

    for PPI in PPIs:
        proteinA = PPI[0]
        proteinB = PPI[1]

        # remove protein hubs
        if (proteinA in proteinHubs) or (proteinB in proteinHubs):
            continue
        else:
            if (proteinA in Uniprot2ENSG) and (proteinB in Uniprot2ENSG):
                # remove self-loops
                if proteinA == proteinB:
                    continue
                else:
                    out_line = (Uniprot2ENSG[proteinA], "pp", Uniprot2ENSG[proteinB])
                    print('\t'.join(out_line))

    logging.info("All done, completed successfully!")

    return


def main():
    file_parser = argparse.ArgumentParser(description="""
                                          Parses the output file(s) produced by interaction_parser.py
                                          and uniprot_parsed.tsv produced by uniprot_parser.py
                                          to produce a high-quality interactome (ie, no hubs and no self-loops)
                                          in SIF format:
                                          - ENSG of Protein A
                                          - "pp" for "protein-protein interaction"
                                          - ENSG of Protein B
                                          """,
                                          formatter_class=argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')

    required.add_argument('--interactions_parsed_files',
                          nargs='+',
                          required=True)
    required.add_argument('--uniprot_file',
                          required=True)

    args = file_parser.parse_args()
    Interactome_Uniprot2ENSG(args)


if __name__ == "__main__":
    # Logging to Standard Error
    logging.basicConfig(format="%(levelname)s %(asctime)s: %(filename)s - %(message)s",
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        level=logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I')
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
