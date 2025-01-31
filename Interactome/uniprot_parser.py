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
    - Gene Name (or a comma seperated list)
    """

    logging.info("Starting to run")

    uniprot_file = sys.stdin

    header = ['Primary_AC', 'TaxID', 'ENSTs', 'ENSGs', 'Secondary_ACs', 'GeneIDs', 'GeneNames']
    print('\t'.join(header))

    ACs = ''
    taxID = 0
    ENSTs = []
    ENSGs = []
    geneIDs = []
    geneNames = []

    # Accession Numbers (AC), strip trailing ';'
    re_AC = re.compile(r'^AC\s+(\S.*);$')

    # get Gene Name / Synonyms from GN line,
    # some records are incomplete and do not end with ";"
    # or unclosed curly braces
    # ex: GN   Name=aadK {ECO:0000303|PubMed:17609790, ECO:0000303|PubMed:8293959,
    # sometimes, GN lines starts with 'Synonyms'
    # ex: GN   Synonyms=AO {ECO:0000303|PubMed:27255930}
    re_GN = re.compile(r'^GN\s+((Name=\S.*)|(Synonyms=\S.*))[,;]$')

    # get TaxID from OX line; some lines violate the Uniprot spec
    # ex: NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
    re_taxID = re.compile(r'^OX\s+NCBI_TaxID=(\d+)[; ]')

    # get Ensembl from the DR line
    re_Ensembl = re.compile(r'^DR\s+Ensembl; (\w+).*; \S+; (\w+).*\.')

    # get GeneIDs from the DR line
    re_geneID = re.compile(r'^DR\s+GeneID;\s+(\d+);')

    for line in uniprot_file:
        line = line.rstrip('\r\n')

        if (re_AC.match(line)):
            # separate ACs
            if (ACs != ''):
                ACs += '; '
            ACs += re_AC.match(line).group(1)
        elif (re.match(r'^AC\s', line)):
            logging.error("Missed the AC line: \n" + line)
            sys.exit()
        elif (re_GN.match(line)):
            GN_line = re_GN.match(line).group(1)
            try:
                # it GN line contains info other than Gene Name,
                # it will also be seperated by a '; '
                # ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
                GN_line_split = GN_line.split('; ')
                if GN_line_split:
                    for i in range(len(GN_line_split)):
                        if re.match(r'^Name=(\S.*)', GN_line_split[i]):
                            geneInfo = re.match(r'^Name=(\S.*)', GN_line_split[i]).group(1)
                            # remove additional info (eg. Pubmed ID) from Gene Name
                            # ex: GN   Name=dbaA {ECO:0000303|PubMed:23001671};
                            try:
                                # split at ' {' and keep only Gene Name
                                geneName = geneInfo.split(' {')[0]
                            except Exception:
                                geneName = geneInfo

                            if geneName not in geneNames:
                                geneNames.append(geneName)

                        # get Synonyms for Gene Name
                        if re.match(r'^Synonyms=(\S.*)', GN_line_split[i]):
                            geneSynonymsInfo = re.match(r'^Synonyms=(\S.*)', GN_line_split[i]).group(1)
                            geneSynonyms = []
                            # there can be multiple synonyms seperated by a ','
                            # ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da;
                            try:
                                geneSynonymList = geneSynonymsInfo.split(', ')
                                if geneSynonymList:
                                    # remove additional info from Synonyms
                                    # ex: GN   Name=Sh3bp5; Synonyms=Sab {ECO:0000303|PubMed:10339589};
                                    for synonym in geneSynonymList:
                                        try:
                                            geneSynonym = synonym.split(' {')[0]
                                            geneSynonyms.append(geneSynonym)
                                        except Exception:
                                            geneSynonyms.append(synonym)
                            except Exception:
                                geneSynonyms.append(geneSynonymsInfo)
                            for synonym in geneSynonyms:
                                # remove additional info
                                if ':' not in synonym:
                                    if synonym not in geneNames:
                                        geneNames.append(synonym)
            except Exception:
                # if GN line contains only Gene Name, do not split
                # ex: GN   Name=APOM;
                if re.match(r'^Name=(\S.*)', GN_line):
                    geneInfo = re.match(r'^Name=(\S+)', GN_line).group(1)
                    try:
                        geneNamewithAcID = geneInfo.split(' {')
                        geneName = geneNamewithAcID[0]
                    except Exception:
                        geneName = geneInfo
                    if geneName not in geneNames:
                        geneNames.append(geneName)
        elif (re.match(r'^GN\s+Name.*', line)):
            logging.error("Missing Gene Name \n" + ACs + line)
            sys.exit()
        elif (re.match(r'^GN\s+Synonyms.*', line)):
            logging.error("Missing Gene Name Synonym \n" + ACs + line)
            sys.exit()
        elif (re_taxID.match(line)):
            if (taxID != 0):
                logging.error("Several OX lines for protein: \t" + ACs)
                sys.exit()
            taxID = re_taxID.match(line).group(1)
        elif (re.match(r'^OX\s', line)):
            logging.error("Missing OX line \n" + line)
            sys.exit()
        elif (re_Ensembl.match(line)):
            ENST_match = re_Ensembl.match(line)

            ENST = ENST_match.group(1)
            if ENST not in ENSTs:
                ENSTs.append(ENST)

            ENSG = ENST_match.group(2)
            if ENSG not in ENSGs:
                ENSGs.append(ENSG)

        elif (re.match(r'^DR\s+Ensembl;', line)):
            logging.error("Failed to get all the Ensembl IDs \n" + ACs + line)
            sys.exit()
        elif (re_geneID.match(line)):
            geneIDs.append(re_geneID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)):
            logging.error("Missing GeneIDs \n" + ACs + line)
            sys.exit()

        elif (line == '//'):  # '//' == end of the record
            if (taxID == '9606'):  # 9606 == human
                try:
                    ACs_split = ACs.split('; ')
                    primAC = ACs_split[0]
                    secACs = ','.join(ACs_split[1:])
                except Exception:
                    logging.error("Failed to store Uniprot Accession IDs for protein: \t" + ACs)
                    sys.exit()

                try:
                    geneNames = ','.join(geneNames)
                except Exception:
                    logging.error('Failed to store Gene Name for protein: \t' + ACs)
                    sys.exit()

                try:
                    ENSTs = ','.join(ENSTs)
                    ENSGs = ','.join(ENSGs)
                except Exception:
                    logging.error("Failed to store ENST/ENSG for protein: \t" + ACs)
                    sys.exit()

                try:
                    geneIDs = ','.join(geneIDs)
                except Exception:
                    raise Exception("Failed to store GeneID for protein: \t' + ACs")

                # print to STDOUT
                out_line = [primAC, taxID, ENSTs, ENSGs, secACs, geneIDs, geneNames]
                print('\t'.join(out_line))

            # reset and move on to the next record
            ACs = ''
            taxID = 0
            ENSTs = []
            ENSGs = []
            geneIDs = []
            geneNames = []

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
