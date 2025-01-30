# This script was cloned from git@github.com:manojmw/grexome-TIMC-Secondary.git

import sys
import re
import argparse
import logging


def parse_uniprot_file(uniprot_file):
    """
    Parse on STDIN a Uniprot file and extracts the following fields from each record:
    - Primary Accession and Secondary Accession from the 'AC' line
    - Gene name and synonyms from the 'GN' line
    - Taxonomy identifier from the 'OX' line
    - ENST(s), ENSG(s) & GeneID(s) from the 'DR' line

    Print to STDOUT in TSV format:
    - Uniprot Primary Accession
    - Taxonomy identifier
    - ENST (or a comma seperated list)
    - ENSG (or a comma seperated list)
    - Uniprot Secondary Accession (or a comma seperated list)
    - GeneID (or a comma seperated list)
    - Gene Name (or a comma seperated list of)
    """

    logging.info("Starting to run")

    uniprot_file = sys.stdin

    header = ['Primary_AC', 'TaxID', 'ENSTs', 'ENSGs', 'Secondary_ACs', 'GeneIDs', 'GeneNames']
    print('\t'.join(header))

    ACs = ''
    TaxID = 0
    ENSTs = []
    ENSGs = []
    GeneIDs = []
    GeneNames = []

    # Accession Numbers (AC), strip trailing ';'
    re_AC = re.compile(r'^AC\s+(\S.*);$')

    # Gene Name/Synonyms from the GN line
    # the below regex pattern ends with either ; or , because
    # UniProt has a bug (at the time of writing this script)
    # the lines for some records are incomplete and do not
    # end with ; as supposed to and even unclosed curly braces
    # Ex: GN   Name=aadK {ECO:0000303|PubMed:17609790, ECO:0000303|PubMed:8293959,
    # Sometimes, GN lines can also start just with 'Synonyms'
    # Ex: GN   Synonyms=AO {ECO:0000303|PubMed:27255930}
    re_GN = re.compile(r'^GN\s+((Name=\S.*)|(Synonyms=\S.*))[,;]$')

    # TaxID from the OX line; some lines violate the Uniprot spec
    # Grab only TaxID even if additional info comes after the TaxID
    # Ex: NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
    re_TaxID = re.compile(r'^OX\s+NCBI_TaxID=(\d+)[; ]')

    # Ensembl transcripts and genes from the DR line
    re_ENS = re.compile(r'^DR\s+Ensembl; (\w+).*; \S+; (\w+).*\.')

    # GeneIDs from the DR line
    re_GID = re.compile(r'^DR\s+GeneID;\s+(\d+);')

    for line in uniprot_file:
        line = line.rstrip('\r\n')

        if (re_AC.match(line)):
            if (ACs != ''):
                # Add trailing separator to the previous entries - distingushes each AC
                ACs += '; '
            ACs += re_AC.match(line).group(1)
        elif (re.match(r'^AC\s', line)):
            # If any AC line is missed, Exit the program with an error message
            logging.error("Missed the AC line: \n" + line)
            sys.exit()
        elif (re_GN.match(line)):
            GNLine = re_GN.match(line).group(1)
            # Grab everything from the GN line
            # The GN line always ends with a ';' except in some cases
            # Ex: GN   Name=aadK {ECO:0000303|PubMed:17609790, ECO:0000303|PubMed:8293959,
            try:
                # If the GN line contains any additional info other than Gene Name,
                # it will also be seperated by a '; '
                # Ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
                GNList = GNLine.split('; ')
                if GNList:
                    for geneinfo in GNList:
                        # Retrieving only the Gene Name
                        if re.match(r'^Name=(\S.*)', geneinfo):
                            GeneNameD = re.match(r'^Name=(\S.*)', geneinfo).group(1)
                            # Remove additional info (eg, Pubmed id) from Gene Name
                            # Ex: GN   Name=dbaA {ECO:0000303|PubMed:23001671};
                            try:
                                # So we split at the ' {' and keep only Gene Name
                                GeneNamewithAcID = GeneNameD.split(' {')
                                GeneName = GeneNamewithAcID[0]
                            except Exception:
                                GeneName = GeneNameD
                                # There can be multiple 'GN' lines, especially
                                # when there are many synonyms
                                # So we check to avoid adding the same gene name again
                                # Ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
                                # Ex: GN   Name=Jon99Ciii; Synonyms=SER2, SER5, Ser99Db; ORFNames=CG15519;
                            if GeneName not in GeneNames:
                                GeneNames.append(GeneName)
                        # retreiving synonyms for Gene Name
                        if re.match(r'^Synonyms=(\S.*)', geneinfo):
                            GeneSynoinfo = re.match(r'^Synonyms=(\S.*)', geneinfo).group(1)
                            GeneSynonyms = []
                            # There can be multiple synonyms seperated by a ','
                            # Ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da;
                            try:
                                GeneSynonymList = GeneSynoinfo.split(', ')
                                if GeneSynonymList:
                                    # Remove additional info from Gene Synonyms
                                    # Ex: GN   Name=Sh3bp5; Synonyms=Sab {ECO:0000303|PubMed:10339589};
                                    for synonym in GeneSynonymList:
                                        try:
                                            Gsynowithaddinfo = synonym.split(' {')
                                            Gsynonym = Gsynowithaddinfo[0]
                                            GeneSynonyms.append(Gsynonym)
                                        except Exception:
                                            GeneSynonyms.append(synonym)
                            except Exception:
                                GeneSynonyms.append(GeneSynoinfo)
                            for synonym in GeneSynonyms:
                                # Remove additional info (i.e not a synonym)
                                if ':' not in synonym:
                                    if synonym not in GeneNames:
                                        GeneNames.append(synonym)
            except Exception:
                # if GN line contains only Gene Name, do not split
                # Ex: GN   Name=APOM;
                if re.match(r'^Name=(\S.*)', GNLine):
                    GeneNameD = re.match(r'^Name=(\S+)', GNLine).group(1)
                    try:
                        # When Gene Name contains additional info
                        GeneNamewithAcID = GeneNameD.split(' {')
                        GeneName = GeneNamewithAcID[0]
                    except Exception:
                        GeneName = GeneNameD
                    if GeneName not in GeneNames:
                        GeneNames.append(GeneName)
        elif (re.match(r'^GN\s+Name.*', line)):
            logging.error("Error: Missed the Gene Name \n" + ACs + line)
            sys.exit()
        elif (re.match(r'^GN\s+Synonyms.*', line)):
            logging.error("Error: Missed the Gene Name Synonym \n" + ACs + line)
            sys.exit()
        elif (re_TaxID.match(line)):
            if (TaxID != 0):
                logging.error("Several OX lines for the protein: \t" + ACs)
                sys.exit()
            TaxID = re_TaxID.match(line).group(1)
        elif (re.match(r'^OX\s', line)):
            logging.error("Missed the OX line \n" + line)
            sys.exit()
        elif (re_ENS.match(line)):
            ENS_match = re_ENS.match(line)
            ENST = ENS_match.group(1)
            if ENST not in ENSTs:
                ENSTs.append(ENST)
            ENSG = ENS_match.group(2)
            if ENSG not in ENSGs:
                ENSGs.append(ENSG)
        elif (re.match(r'^DR\s+Ensembl;', line)):
            logging.error("Failed to get all the Ensembl identifiers \n" + ACs + line)
            sys.exit()
        elif (re_GID.match(line)):
            GeneIDs.append(re_GID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)):
            logging.error("Missed the GeneIDs \n" + ACs + line)
            sys.exit()

        elif (line == '//'):  # '//' end of the record
            if (TaxID == '9606'):  # Human
                try:
                    ACs_split = ACs.split('; ')
                    primary_AC = ACs_split[0]
                    secondary_AC_list = ACs_split[1:]
                    secondary_ACs = ','.join(secondary_AC_list)
                except Exception:
                    logging.error("Failed to store store UniProt Accession IDs for the protein: \t" + ACs)
                    sys.exit()

                try:
                    GeneNames = ','.join(GeneNames)
                except Exception:
                    logging.error('Error: Failed to store Gene Name for the protein: \t' + ACs)
                    sys.exit()

                try:
                    ENSTs = ','.join(ENSTs)
                    ENSGs = ','.join(ENSGs)
                except Exception:
                    logging.error("Failed to store ENST / ENSG for the protein: \t" + ACs)
                    sys.exit()

                try:
                    GeneIDs = ','.join(GeneIDs)
                except Exception:
                    raise Exception("Error: Failed to store Gene identifiers for the protein: \t' + ACs")

                # Print to STDOUT
                UniProt_outline = [primary_AC, TaxID, ENSTs, ENSGs, secondary_ACs, GeneIDs, GeneNames]
                print('\t'.join(UniProt_outline))

            # Reset and move on to the next record
            ACs = ''
            TaxID = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []
            GeneNames = []

            continue

    logging.info("All done, completed successfully!")

    return


def main():
    file_parser = argparse.ArgumentParser(description="""
                                          Parse on STDIN a Uniprot file and print to STDOUT in TSV format:
                                          - Uniprot Primary Accession
                                          - Taxonomy identifier
                                          - ENST (or a comma seperated list)
                                          - ENSG (or a comma seperated list)
                                          - Uniprot Secondary Accession (or a comma seperated list)
                                          - GeneID (or a comma seperated list)
                                          - Gene Name (or a comma seperated list)
                                          """,
                                          formatter_class=argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=parse_uniprot_file)
    args = file_parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s %(asctime)s: %(filename)s - %(message)s",
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        level=logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I')
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
