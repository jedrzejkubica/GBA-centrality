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
import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_uniprot_file(uniprot_file):
    """
    Parse on STDIN a Uniprot file.

    Print to STDOUT in TSV format:
    - Uniprot Primary AC
    - Uniprot Secondary AC(s)
    - tax ID
    - gene name(s)
    - ENSG(s)
    """

    # print header to STDOUT
    print("\t".join(["PrimaryAC", "SecondaryACs", "TaxID", "GeneNames", "ENSGs"]))

    re_AC = re.compile(r'^AC\s+(\S.*);$')  # optional trailing ";"
    # GN lines might have optionally: multiple names separated by ",", "{}" block and trailing ";"
    re_GN = re.compile(r'^GN\s+((Name=\S.*)|(Synonyms=\S.*))[,;]?$')
    re_Name = re.compile(r'^Name=([\S][^{;]+?)(?:\s*\{[^}]*\})?;?$')
    re_Synonyms = re.compile(r'^Synonyms=([\S][^{;]+?)(?:\s*\{[^}]*\})?;?$')

    re_taxID = re.compile(r'^OX\s+NCBI_TaxID=(\d+);?$')
    # ENSGs and ENSTs should be in the same line
    re_Ensembl = re.compile(r'^DR\s+Ensembl;')
    re_ENSG = re.compile(r'(ENSG[\d\.]+)')

    primary_AC = ""
    secondary_ACs = []
    taxID = ""
    gene_names = []
    ENSGs = []

    for line in uniprot_file:
        if re_AC.match(line):
            ACs_split = re_AC.match(line).group(1).split("; ")
            primary_AC = ACs_split[0]
            secondary_ACs = ACs_split[1:]
        elif re_GN.match(line):
            GN_split = re_GN.match(line).group(1).split("; ")
            for el in GN_split:
                if re_Name.match(el):
                    names_split = re_Name.match(el).group(1).split(', ')
                    for name in names_split:
                        gene_names.append(name)
                elif re_Synonyms.match(el):
                    synonyms_split = re_Synonyms.match(el).group(1).split(', ')
                    for synonym in synonyms_split:
                        gene_names.append(synonym)
        elif re_taxID.match(line):
            taxID = re_taxID.match(line).group(1)
        elif re_Ensembl.match(line):
            Ensembl_split = line.split("; ")
            for el in Ensembl_split:
                if re_ENSG.match(el):
                    ENSG = re_ENSG.match(el).group(1).split('.')[0]  # remove version numbers
                    ENSGs.append(ENSG)
        elif line.startswith("//"):  # end of the record
            if(primary_AC != "" and taxID != "" and len(gene_names) != 0 and len(ENSGs) != 0):
                out_line = [primary_AC,
                            ",".join(secondary_ACs),
                            taxID,
                            ",".join(gene_names),
                            ",".join(ENSGs)]
                print('\t'.join(out_line))
                
            primary_AC = ""
            secondary_ACs = []
            taxID = ""
            gene_names = []
            ENSGs = []

            continue


def main():

    logger.info("Parsing uniprot file")
    parse_uniprot_file(sys.stdin)

    logger.info("Done!")

if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logger, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    try:
        main()
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
