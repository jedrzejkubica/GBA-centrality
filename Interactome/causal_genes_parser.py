import os
import sys
import argparse
import logging
import re

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_uniprot(uniprot_file):
    '''
    Parses tab-seperated Uniprot file produced by Interactome/uniprot_parser.py
    which consists of 7 columns (one record per line):
    - Uniprot Primary Accession
    - Taxonomy Identifier
    - ENST (or a comma seperated list of ENSTs)
    - ENSG (or a comma seperated list of ENSGs)
    - Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
    - GeneID (or a comma seperated list of GeneIDs)
    - Gene Name (or a comma seperated list of Gene Names)

    Returns:
      - gene2ENSG: dict with key=gene, value=ENSG
      - ENSG2uniprot: dict with key=ENSG, value=Primary accession

    Note: if more than one gene name is associated with a particular ENSG,
          then keeping the first gene name from the list
    '''
    gene2ENSG = {}
    ENSG2uniprot = {}

    try:
        f = open(uniprot_file, 'r')
    except Exception as e:
        logging.error("Opening provided uniprot file %s: %s", uniprot_file, e)
        raise Exception("cannot open provided Uniprot file")

    # skip header
    line = f.readline()
    if not line.startswith("Primary_AC\t"):
        logging.error("uniprot file %s is headerless? expecting headers but got %s",
                      uniprot_file, line)
        raise Exception("Uniprot file problem")

    for line in f:
        split_line = line.rstrip('\r\n').split('\t')

        # if some records are incomplete, die
        if len(split_line) != 7:
            logging.error("uniprot file %s line doesn't have 7 fields: %s",
                          uniprot_file, line)
            raise Exception("Uniprot file problem")

        (AC_primary, TaxID, ENSTs, ENSGs, AC_secondary, GeneIDs, geneNames) = split_line

        # make sure there is at least one ENSG and keep only the first one
        if ENSGs == "":
            continue
        ENSG = ENSGs.split(',')[0]

        # make sure there is at least one gene name and keep only the first one
        if geneNames == "":
            continue
        geneName = geneNames.split(',')[0]

        gene2ENSG[geneName] = ENSG
        ENSG2uniprot[ENSG] = AC_primary

    return(gene2ENSG, ENSG2uniprot)


def parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2uniprot):
    '''
    Build a list of protein Uniprot Primary AC corresponding to
    causal gene names from causal_genes_file

    arguments:
    - causal_genes_file: filename (with path) of known causal genes, one gene name per line
    - gene2ENSG: dict of all known genes, key=gene_name, value=ENSG
    - ENSG2uniprot: type=dict, key=ENSG, value=Primary accession

    returns:
    - causal_proteins: list of Uniprot Primary accession for causal genes
    '''
    causal_proteins = []
    num_found_genes = 0

    try:
        f_causal = open(causal_genes_file, 'r')
    except Exception as e:
        logger.error("Opening provided causal genes file %s: %s", causal_genes_file, e)
        raise Exception("cannot open provided causal genes file")

    re_causal = re.compile(r'^[a-zA-Z0-9\-_]+$')  # allow for: letters, digits, "_", "-"

    for line in f_causal:
        if re_causal.match(line):
            gene_name = line.rstrip()
            if gene_name in gene2ENSG:
                ENSG = gene2ENSG[gene_name]
                if ENSG in ENSG2uniprot:
                    causal_proteins.append(ENSG2uniprot[ENSG])
                    num_found_genes += 1
                else:
                    logger.warning("causal gene %s == %s is not in Uniprot, skipping it",
                                   gene_name, ENSG)
            else:
                logger.warning("causal gene %s is not a known gene in gene2ENSG, skipping it",
                               gene_name)
        else:
            logger.error("Bad line in the causal genes file, doesn't look like a gene name: %s", line)
            raise Exception("Bad line in the causal genes file")

    f_causal.close()

    logger.info("found %i causal genes with known ENSG", num_found_genes)

    return(causal_proteins)


def save_causal(causal_proteins):
    '''
    Print scores to stdout in TSV format, 2 columns: node score
    a seeds.txt file with one protein per line

    arguments:
    - causal_proteins: list of causal proteins
    '''
    for protein in causal_proteins:
        print(protein)


def main(uniprot_file, causal_genes_file):

    logger.info("Parsing Uniprot file")
    (gene2ENSG, ENSG2uniprot) = parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    (causal_proteins) = parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2uniprot)

    logger.info("Printing protein seeds")
    save_causal(causal_proteins)

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
        Parses the causal genes file and uniprot_parsed.tsv produced by uniprot_parser.py,
        print causal proteins to stdout
        """)
    parser.add_argument('--uniprot', required=True)
    parser.add_argument('--causal', required=True)

    args = parser.parse_args()

    try:
        main(args.uniprot, args.causal)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)