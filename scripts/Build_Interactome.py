# This script was cloned from git@github.com:manojmw/grexome-TIMC-Secondary.git

#!/usr/bin/python

# manojmw
# 27 Jan, 2022

import sys, argparse
import logging
import gzip

###########################################################

# Parses the tab-seperated output files
# produced by the Interaction_parser.py (No headers)
#
# Required columns are:
# - Protein A UniProt PrimAC (column-1)
# - Protein B UniProt PrimAC (column-2)
# - Interaction Detection Method (column-3)
# - Pubmed Identifier (column-4)
# - Interaction type (column-5)
#
# Processes it by filtering based on Interaction Detection Method and Interaction type
# Final filtering - Each interaction has at least 2 experiments of which
# at least one should be proven by any binary interaction detection method
#
# Returns a list with 5 items (in each sublist):
# - Protein A UniProt PrimAC,
# - Protein B UniProt PrimAC,
# - Publication_Count,
# - PMID(s) and
# - Experiment_count
def UniProtInteractome(inExpFile):

    logging.info("starting to run")

    # Dictionary for storing Interacting
    # Proteins and Pubmed Identifier
    PPI_PMID_dict = {}

    # Dictionary for storing Interacting
    # Proteins and Interaction Detection Method
    PPI_IntDetMethod_dict = {}

     # List of User input PPI interaction experiments file(s)
    PPIExpFiles = inExpFile

    # there can be multiple files
    for file in PPIExpFiles:

        PPIExpFile = open(file)

        # Data lines
        for line in PPIExpFile:

            line = line.rstrip('\n')

            ExpPPI_fields = line.split('\t')

            # Filtering out Interactions based on Interaction Detection Method
            IntDetMethod = ExpPPI_fields[2]
            # MI:0096 -> pull down
            # MI:0254 -> genetic interference
            # MI:0686 -> unspecified method

            ###Filtering out Interactions based on Interaction Type
            IntType = ExpPPI_fields[4].rstrip('\n')
            # MI:0407 -> direct interaction
            # MI:0915 -> physical association

            if IntDetMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and IntType in ['MI:0407', 'MI:0915']:
                # ExpPPI_fields[0] -> Protein_A_UniprotPrimAC
                # ExpPPI_fields[1] -> Protein_B_UniprotPrimAC
                Interactors = ExpPPI_fields[0] + '_' + ExpPPI_fields[1]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Pubmed Identifier (PMID) - ExpPPI_fields[3]
                (Int_key, PMIDs) = (Interactors, ExpPPI_fields[3])

                # Check if the Key exists in PPI_PMID_dict
                # If yes, then store the values (PMIDs) as a list
                if PPI_PMID_dict.get(Int_key, False):
                    # Avoiding duplicate PMIDs
                    if PMIDs not in PPI_PMID_dict[Int_key]:
                        PPI_PMID_dict[Int_key].append(PMIDs)
                else:
                    PPI_PMID_dict[Int_key] = [PMIDs]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Interaction Detection Method - ExpPPI_fields[2]
                (Int_key, IntDetMeth) = (Interactors, ExpPPI_fields[2])

                if PPI_IntDetMethod_dict.get(Int_key, False):
                    PPI_IntDetMethod_dict[Int_key].append(IntDetMeth)
                else:
                    PPI_IntDetMethod_dict[Int_key] = [IntDetMeth]

        # Closing the file
        PPIExpFile.close()

    # Initializing output list
    Uniprot_Interactome_list = []

    # Processing the dictionaries and returning a list
    # Since both the dictionaries have the same keys,
    # We can use the keys from PPI_PMID_dict to iterate through PPI_IntDetMethod_dict
    for Int_key in PPI_PMID_dict:

        # Checking if at least one of the experiments for a given interaction has been proved by any binary interaction method
        # i.e. Other than Affintity Chromatography Technology (ACT) - 'MI:0004'
        # Used to eliminate PPIs that have been proved ONLY by using ACT
        if any(value != 'MI:0004' for value in PPI_IntDetMethod_dict[Int_key]):
            Proteins = Int_key.split('_')
            Protein_A = Proteins[0]
            Protein_B = Proteins[1]

            Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key])
            PMID_count = len(PPI_PMID_dict[Int_key])
            Exp_count = len(PPI_IntDetMethod_dict[Int_key])

            if Exp_count >= 2:
                interaction_out_line = [Protein_A, Protein_B, str(PMID_count), Pubmed_Identifier, str(Exp_count)]
                Uniprot_Interactome_list.append(interaction_out_line)

    return Uniprot_Interactome_list

###########################################################

# Parses tab-seperated Uniprot file produced by Uniprot_parser.py
# which consists of 7 columns (one record per line):
# - Uniprot Primary Accession
# - Taxonomy Identifier
# - ENST (or a comma seperated list of ENSTs)
# - ENSG (or a comma seperated list of ENSGs)
# - Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
# - GeneID (or a comma seperated list of GeneIDs)
# - Gene Name (or a comma seperated list of Gene Names)
#
# Returns:
#   - ENSG_Gene_dict: dict with key=ENSG, values=geneName
#   - Uniprot_ENSG_dict: dict with key=Primary accession, values=ENSG

# Note: if more than one gene name is associated with a particular ENSG,
#       then keeping the first gene name from the list
def ParseUniprot(inUniProt):
    ENSG_Gene_dict = {}
    Uniprot_ENSG_dict = {}

    try:
        f = open(inUniProt, 'r')
    except Exception as e:
        logging.error("Failed to read the Uniprot file: %s" % inUniProt, e)
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

        ENSG_Gene_dict[ENSG] = geneName
        Uniprot_ENSG_dict[AC_primary] = ENSG

    return ENSG_Gene_dict, Uniprot_ENSG_dict

###########################################################

# Parses the Uniprot_Interactome_list returned by 
# the function UniProtInteractome
#
# Required items in each sublist:
# UniProt Primary Accession of Protein A (first item)
# UniProt Primary Accession of Protein B (second item)
#
# Returns 2 dictionaries
# First dictionary contains:
# - key: Protein A; Value: List of interactors
# Second dictionary contains:
# - key: Protein B; Value: List of interactors
def Interactome_dict(Uniprot_Interactome_list):

    # Dictionaries to store interacting proteins
    # In ProtA_dict, key -> Protein A; Value -> Protein B
    # In ProtB_dict, key -> Protein B; Value -> Protein A
    ProtA_dict = {}
    ProtB_dict = {}
        
    # Data lines
    for data in Uniprot_Interactome_list:

        if data[0] != data[1]: # checking self-interaction
            # Check if the Key(ProtA) exists in ProtA_dict
            # If yes, then append the interctor to 
            # the list of values (Interactors)
            if ProtA_dict.get(data[0], False):
                ProtA_dict[data[0]].append(data[1])
            else:
                ProtA_dict[data[0]] = [data[1]]

            # Check if the Key(ProtB) exists in ProtB_dict
            # If yes, then append the interctor to 
            # the list of values (Interactors)
            if ProtB_dict.get(data[1], False):
                ProtB_dict[data[1]].append(data[0])
            else:
                ProtB_dict[data[1]] = [data[0]]    

        # else:
            # NOOP -> The interaction is a self-interaction

    return ProtA_dict, ProtB_dict

###########################################################

# Parses the ProtA_dict & ProtB_dict dictionaries
# returned by the function: Interactome_dict
#
# Checks if a protein is hub/sticky protein
# The criteria for considering a given protein
# as a hub is if it interacts with > 120 proteins
# This number is based on the degree distribution
# of all the proteins in the entire high-quality Interactome
# before eliminating hub/sticky proteins
#
# Returns a dictionary:
# - Key: UniProt Primary Accession of the Hub protein
# - Value: 1
def getHubProteins(ProtA_dict, ProtB_dict):

    # Initializing dictionary to store hub proteins
    HubProteins = {}

    # Checking the number of Interactors for a 
    # given protein
    for protein in ProtA_dict:
        # Get the no. of Interactors
        InteractorsCount_ProtA_dict = len(ProtA_dict[protein])
        # Check if this protein is also
        # present in ProtB_dict
        if protein in ProtB_dict:
            InteractorsCount_ProtB_dict = 0
            # If present, loop through each interactor to 
            # make sure that the Interactors
            # for the current protein in ProtB_dict was not already seen
            # in the Interactor list of ProtA_dict
            for interactor in ProtB_dict[protein]:
                if not interactor in ProtA_dict[protein]:
                    InteractorsCount_ProtB_dict += 1
            Total_InteractorsCount = InteractorsCount_ProtA_dict + InteractorsCount_ProtB_dict
        # if the protien not present in ProtB_dict
        # simply get the Interactors count from ProtA_dict
        else:
            Total_InteractorsCount = InteractorsCount_ProtA_dict

        # If the protein has > 120 Interactors 
        # it is considered a hub/sticky protein
        # append it to the HubProteins list
        if Total_InteractorsCount > 120:
            HubProteins[protein] = 1

    # Checking the interactor count of proteins
    # present only in ProtB_dict
    for protein in ProtB_dict:
        if not protein in ProtA_dict:
            # checking for Hub protein
            if len(ProtB_dict[protein]) > 120:
                HubProteins[protein] = 1

    return HubProteins 

###########################################################

# Parses the Uniprot_Interactome_list
# returned by the function: UniProtInteractome()
# and the dictionaries returned by ParseUniprot()
#
# Prints to STDOUT in SIF format
# Output consists of 2 columns:
# - ENSG of Protein A
# - "pp" for "protein-protein interaction"
# - ENSG of Protein B
#
# Note: Hub proteins and self-loops are removed
def Interactome_Uniprot2ENSG(args):

    # Calling the functions
    Uniprot_Interactome_list = UniProtInteractome(args.inExpFile)
    ENSG_Gene_dict, Uniprot_ENSG_dict = ParseUniprot(args.inUniProt)
    (ProtA_dict, ProtB_dict) = Interactome_dict(Uniprot_Interactome_list)
    HubProteins = getHubProteins(ProtA_dict, ProtB_dict)

    for sublist in Uniprot_Interactome_list:
        proteinA = sublist[0]
        proteinB = sublist[1]

        # Removing hub/sticky proteins
        if (proteinA in HubProteins) or (proteinB in HubProteins):
            continue 
        else:

            if (proteinA in Uniprot_ENSG_dict) and (proteinB in Uniprot_ENSG_dict):
                # Removing self-loops
                if proteinA == proteinB:
                    continue
                else:
                    ENSG_Interactome_out = (Uniprot_ENSG_dict[proteinA], "pp", Uniprot_ENSG_dict[proteinB])
                    print('\t'.join(ENSG_Interactome_out))

    logging.info("ALL DONE, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
    --------------------------------------------------------------------------------------------------
    Parses the output file(s) produced by Interaction_parser.py 
    and Uniprot_output.tsv produced by Uniprot_parser.py
    to produce a high-quality interactome (ie, no hubs, no self-loops).
    --------------------------------------------------------------------------------------------------
    The output (high-quality interactome) consists of three columns in SIF format:
    - ENSG of Protein A
    - "pp" for "protein-protein interaction"
    - ENSG of Protein B
    --------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')

    required.add_argument('--inExpFile', metavar = "Input File", dest = "inExpFile", nargs = '+', help = 'Output files produced by Interaction_parser.py', required = True)
    required.add_argument('--inUniProt', metavar = "Input File", dest = "inUniProt", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)

    args = file_parser.parse_args()
    Interactome_Uniprot2ENSG(args)


if __name__ == "__main__":
    # Logging to Standard Error
    logging.basicConfig(format = "%(levelname)s %(asctime)s: %(filename)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S', stream = sys.stderr, level = logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I' )
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
