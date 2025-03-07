# This script was cloned from git@github.com:manojmw/grexome-TIMC-Secondary.git

import os
import sys
import re
import argparse
import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_uniprot_file(uniprot_file):
    """
    Parse the TSV file from uniprot_parser.py with columns:
    - Primary_AC
    - TaxID
    - Secondary_ACs
    - GeneIDs
    - GeneNames

    Returns:
    - primAC2taxID: dict with key=Uniprot Primary Accession, value=TaxID
    - secAC2primAC: dict with  key=Uniprot Secondary Accession, value=Uniprot Primary Accession
    - geneID2primAC: dict with key=GeneID, value=Uniprot Primary Accession
    - geneName2primAC: dict with key=Gene Name, value=Uniprot Primary Accession
    """

    primAC2taxID = {}
    secAC2primAC = {}
    geneID2primAC = {}
    geneName2primAC = {}

    f = open(uniprot_file)

    header = f.readline()
    header = header.rstrip("\n")
    header_fields = header.split("\t")
    (primAC_idx, taxID_idx, secAC_idx, geneID_idx, geneName_idx) = (-1, -1, -1, -1, -1)

    for i in range(len(header_fields)):
        if header_fields[i] == "Primary_AC":
            primAC_idx = i
        elif header_fields[i] == "TaxID":
            taxID_idx = i
        elif header_fields[i] == "Secondary_ACs":
            secAC_idx = i
        elif header_fields[i] == "GeneIDs":
            geneID_idx = i
        elif header_fields[i] == "GeneNames":
            geneName_idx = i

    if not primAC_idx >= 0:
        logger.error("Missing required column title 'Primary_AC' in the file: \n" + uniprot_file)
        sys.exit()
    elif not taxID_idx >= 0:
        logger.error("Missing required column title 'TaxID' in the file: \n" + uniprot_file)
        sys.exit()
    elif not secAC_idx >= 0:
        logger.error("Missing required column title 'Secondary_ACs' in the file: \n" + uniprot_file)
        sys.exit()
    elif not geneID_idx >= 0:
        logger.error("Missing required column title 'GeneIDs' in the file: \n" + uniprot_file)
        sys.exit()
    elif not geneName_idx >= 0:
        logger.error("Missing required column title 'GeneNames' in the file: \n" + uniprot_file)
        sys.exit()

    for line in f:
        line_split = line.rstrip("\n").split("\t")

        primAC2taxID[line_split[primAC_idx]] = line_split[taxID_idx]

        try:
            secACs = line_split[secAC_idx].split(",")
        except Exception:
            secACs = [line_split[secAC_idx]]

        for i in range(len(secACs)):
            secAC, primAC = (secACs[i], line_split[primAC_idx])

            # if Secondary AC is associated with multiple Primary AC,
            # assign "-1" to avoid using it later
            if secAC2primAC.get(secAC, False):
                if secAC2primAC[secAC] != "-1":
                    secAC2primAC[secAC] = "-1"
            else:
                secAC2primAC[secAC] = primAC

        try:
            geneIDs = line_split[geneID_idx].split(",")
        except Exception:
            geneIDs = [line_split[geneID_idx]]

        for i in range(len(geneIDs)):
            (geneID, primAC) = (geneIDs[i], line_split[primAC_idx])

            # if GeneID is associated with multiple Primary AC,
            # assign "-1" to avoid using it later
            if geneID2primAC.get(geneID, False):
                if geneID2primAC[geneID] != "-1":
                    geneID2primAC[geneID] = "-1"
            else:
                geneID2primAC[geneID] = primAC

        try:
            geneNames = line_split[geneName_idx].split(",")
        except Exception:
            geneNames = [line_split[geneName_idx]]

        for i in range(len(geneNames)):
            (geneName, primAC) = (geneNames[i], line_split[primAC_idx])

            # if GeneName is associated with multiple Primary AC,
            # assign "-1" avoid using it later
            if geneName2primAC.get(geneName, False):
                if geneName2primAC[geneName] != "-1":
                    geneID2primAC[geneName] = "-1"
            else:
                geneName2primAC[geneName] = primAC

    f.close()

    return primAC2taxID, secAC2primAC, geneID2primAC, geneName2primAC


def parse_interaction_file(interaction_file, primAC2taxID, secAC2primAC, geneID2primAC, geneName2primAC):
    """
    Map a miTAB 2.5 or 2.7 file to the parsed Uniprot file from parse_uniprot_file().

    For each interaction, grab Primary Accessions of interacting proteins,
    interaction detection method, PubmedID and interaction type.
    Ignore interactions if:
    - Primary Accessions of either of the proteins cannot be found OR
    - if PubmedID is missing OR
    - if TaxID of both interacting proteins is not "9606" i.e. human

    Print to STDOUT in TSV format:
    - Protein A Uniprot Primary Accession
    - Protein B Uniprot Primary Accession
    - Interaction Detection Method
    - Pubmed ID
    - Interaction Type
    """

    re_uniprot = re.compile(r'^uniprot(kb|/swiss-prot):([A-Z0-9-_]+)$')  # protein Uniprot IDs
    re_uniprot_missed = re.compile(r'^uniprot')

    # if Uniprot AC not found, use GeneID to get the corresponding Primary AC
    re_geneID = re.compile(r'^entrez gene/locuslink:(\d+)$')
    re_geneID_missed = re.compile(r'^entrez gene/locuslink:(\d+)$')

    # if not found using Primary AC, use Gene Name to get the Uniprot Primary AC
    re_geneName = re.compile(r'^(entrez gene.locuslink:|uniprotkb:)([\w\s\_\\\/\:\.\-]+$|[\w\s\_\\\/\:\.\-]+)')

    # grab Interaction Detection Method and Interaction Type
    re_psimi = re.compile(r'^psi-mi:"(MI:\d+)"')
    re_psimi_missed = re.compile(r'^psi-mi:')

    # grab publications that mention the interaction
    re_PubmedID = re.compile(r'^pubmed:(\d+)$')
    re_PubmedID_unassigned = re.compile(r'^pubmed:unassigned')  # some Pubmed IDs are unassigned in IntAct
    re_PubmedID_missed = re.compile(r'^pubmed:')

    interaction_file = open(interaction_file)
    interaction_file.readline()

    for line in interaction_file:
        line_split = line.rstrip("\n").split("\t")

        proteins = ['','']
        detectionMethod = ""
        PubmedID = ""
        interactionType = ""

        # exclude "Complex expansion" column for miTAB 2.7
        # example from IntAct:
        # P81274 & Q14980
        # Interaction Detection Method: 2 hybrid
        # Pubmed:15537540
        # Interaction AC: EBI-624047
        # expansion: spoke expansion
        #
        # exclude "Affinity Chromatography Technology"
        # also corresponds to expansion
        # example from IntAct:
        # P19784 & P03211
        # Interaction Detection Method: Affinity Chromatography Technology
        # Pubmed: 12783858
        # Interaction AC: EBI-8656048
        # expansion: Not spoke expansion
        for protein_idx in [0, 1]:
            if (re_uniprot.match(line_split[protein_idx])):
                ID = re_uniprot.match(line_split[protein_idx]).group(2)

                if primAC2taxID.get(ID, False):
                    proteins[protein_idx] = ID
                    continue
            elif (re_uniprot_missed.match(line_split[protein_idx])):
                logger.error("ID is a Uniprot Accession but failed to grab it for the line:\n", line)
                sys.exit()
            elif (re_geneID.match(line_split[protein_idx])):
                ID = re_geneID.match(line_split[protein_idx]).group(1)
                # if ID == "-1", skip it
                if geneID2primAC.get(ID, "-1") != "-1":
                    proteins[protein_idx] = geneID2primAC[ID]
                    continue
            elif (re_geneID_missed.match(line_split[protein_idx])):
                logger.error("ID is a GeneID but failed to grab it for the line:\n", line)
                sys.exit()

            # if Uniprot AC and GeneID not found,
            # then look in alternative columns in miTAB file;
            # stop at first alternative ID
            altIDs = line_split[2 + protein_idx].split("|")

            for altID in altIDs:
                if (re_uniprot.match(altID)):
                    ID = re_uniprot.match(altID).group(2)

                    if primAC2taxID.get(ID, False):
                        proteins[protein_idx] = ID
                        break
                    # if ID in Secondary ACs, then get Primary AC
                    elif secAC2primAC.get(ID, "-1") != "-1":
                        proteins[protein_idx] = secAC2primAC[ID]
                        break
                elif (re_uniprot_missed.match(altID)):
                    logger.error(altID + " is a Uniprot Accession but failed to grab it for line:\n", line)
                    sys.exit()
                elif (re_geneID.match(altID)):
                    ID = re_geneID.match(altID).group(1)
                    # if ID == "-1", skip it
                    if geneID2primAC.get(ID, "-1") != "-1":
                        proteins[protein_idx] = geneID2primAC[ID]
                        break
                elif (re_geneID_missed.match(altID)):
                    logger.error("alternativeID " + altID + " is a GeneID but failed to grab it for line:\n", line)
                    sys.exit()
                # if Uniprot Primary not found using above, then it using gene name:
                # - in miTAB 2.5, gene name can be in alternative IDs column;
                # - in miTAB 2.7, gene name and gene name synonym can be in the Alias column
                elif (re_geneName.match(altID)):
                    geneNameAltID = re_geneName.match(altID).group(2)
                    # if ID == "-1", skip it
                    if geneName2primAC.get(geneNameAltID, "-1") != "-1":
                        proteins[protein_idx] = geneName2primAC[geneNameAltID]
                        break
            aliasIDs = line_split[4 + protein_idx].split("|")
            for aliasID in aliasIDs:
                if (re_geneName.match(aliasID)):
                    GN_aliasID = re_geneName.match(aliasID).group(2)
                    # if ID == "-1", skip it
                    if geneName2primAC.get(GN_aliasID, "-1") != "-1":
                        proteins[protein_idx] = geneName2primAC[GN_aliasID]
                        break

        # if no Uniprot Primary AC, skip
        if proteins[0] == "" or proteins[1] == "":
            continue

        # grab Interaction Detection Method
        if re_psimi.match(line_split[6]):
            detectionMethod = re_psimi.match(line_split[6]).group(1)
        elif re_psimi_missed.match(line_split[6]):
            logger.error("Failed to grab Interaction Detection Method for line:\n", line)
            sys.exit()

        # grab PubmedID
        # some line include additional fields
        # ex: pubmed:10542231|mint:MINT-5211933, so split at "|"
        PubmedID_fields = line_split[8].split("|")
        for PubmedID_entry in PubmedID_fields:
            if (re_PubmedID.match(PubmedID_entry)):
                PubmedID = re_PubmedID.match(PubmedID_entry).group(1)
            elif (re_PubmedID_unassigned.match(PubmedID_entry)):
                continue
            elif (re_PubmedID_missed.match(PubmedID_entry)):
                logger.error("Failed to grab PubmedID for line:\n", line)
                sys.exit()

        # grab Interaction Type
        if (re_psimi.match(line_split[11])):
            interactionType = re_psimi.match(line_split[11]).group(1)
        elif (re_psimi_missed.match(line_split[11])):
            logger.error("Failed to grab Interaction Type for line:\n", line)
            sys.exit()

        # if no PubmedID, skip
        if (PubmedID == ""):
            continue

        # if not human, skip
        if (primAC2taxID[proteins[0]] != "9606") or (primAC2taxID[proteins[1]] != "9606"):
            continue

        if proteins[0] < proteins[1]:  # sort by Primary AC
            interaction_out_line = [proteins[0], proteins[1], detectionMethod, PubmedID, interactionType]
        else:
            interaction_out_line = [proteins[1], proteins[0], detectionMethod, PubmedID, interactionType]

        print("\t".join(interaction_out_line))

    interaction_file.close()


def main(interaction_file, uniprot_file):

    logger.info("Parsing Uniprot file")
    (primAC2taxID, secAC2primAC, geneID2primAC, geneName2primAC) = parse_uniprot_file(uniprot_file)

    logger.info("Parsing interaction file")
    parse_interaction_file(interaction_file, primAC2taxID, secAC2primAC, geneID2primAC, geneName2primAC)

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
        Parse a miTAB 2.5/2.7 file, then map to Uniprot file
        and print to STDOUT in the TSV format:
        - Uniprot Primary Accession of Protein A
        - Uniprot Primary Accession of Protein B
        - Interaction Detection Method
        - Pubmed Identifier
        - Interaction type
        """)

    parser.add_argument("--interaction_file", required=True)
    parser.add_argument("--uniprot_file", required=True)

    args = parser.parse_args()

    try:
        main(interaction_file=args.interaction_file,
             uniprot_file=args.uniprot_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
