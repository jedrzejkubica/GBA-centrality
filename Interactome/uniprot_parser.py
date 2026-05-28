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
    - Uniprot Secondary AC(s) (comma-separated)
    - tax ID
    - gene name(s) (comma-separated)
    """

    # print header to STDOUT
    print("\t".join(["PrimaryAC", "SecondaryACs", "TaxID", "GeneNames"]))

    re_AC = re.compile(r'^AC\s+(\S.+);$')
    re_taxID = re.compile(r'^OX\s+NCBI_TaxID=([^;]+);$')

    # GN data can be multi-line, parse them all and concatenate, then
    # process at end of record
    re_GN = re.compile(r'^GN\s+(\S.+)$')
    # All gene name(s) and synonym(s) will be output together, comma-separated
    # but without distinguishing between Name, Synonym, or even if multiple
    # genes were listed by uniprot (with "GN  and")

    # REs for extracting Name and Synonyms from concatenated GN data
    re_Name = re.compile(r'Name=([^;]+);')
    re_Synonyms = re.compile(r'Synonyms=([^;]+);')

    # some Names/Synonyms/taxIDs have evidence codes, eg
    # Name=atg-18 {ECO:0000312|WormBase:F41E6.13a};
    # -> re to remove this
    re_removeEC = re.compile(r' {ECO[^}]+}')

    # accumulators for current entry. ACs and taxID get processed on the fly but
    # gene names/synomyms at the end only => store in a string
    primary_AC = ""
    secondary_ACs = []
    taxID = ""
    gene_data = ""

    for line in uniprot_file:
        if re_AC.match(line):
            ACs_split = re_AC.match(line).group(1).split("; ")
            if primary_AC == "":
                primary_AC = ACs_split[0]
                secondary_ACs = ACs_split[1:]
            else:
                secondary_ACs += ACs_split
        elif re_GN.match(line):
            if gene_data != "":
                gene_data += " "
            gene_data += re_GN.match(line).group(1)
        elif re_taxID.match(line):
            taxID = re_taxID.match(line).group(1)
            # remove ECO if present
            taxID = re_removeEC.sub('', taxID)
            # sanity check: should be digits only
            try:
                taxID = int(taxID)
            except Exception:
                raise Exception(f"taxID {taxID} not an integer in {primary_AC}")

        elif line.startswith("//"):  # end of the record
            # process gene_data: extract each name and remove ECO if present
            genes = ""
            names = re_Name.findall(gene_data)
            for name in names:
                if genes != "":
                    genes += ","
                genes += re_removeEC.sub('', name)
            # process symonyms in the same way
            synonymsAll = re_Synonyms.findall(gene_data)
            for synonyms in synonymsAll:
                synonyms_ECremoved = re_removeEC.sub('', synonyms)
                for synonym in synonyms_ECremoved.split(', '):
                    synonym = synonym.strip()
                    genes += ',' + synonym

            if (primary_AC != "" and taxID != ""):
                out_line = [primary_AC,
                            ",".join(secondary_ACs),
                            str(taxID),
                            genes]
                print('\t'.join(out_line))

            primary_AC = ""
            secondary_ACs = []
            taxID = ""
            gene_data = ""

        # sanity check: did we miss anything?
        elif line.startswith(('AC', 'OX', 'GN')):
            logger.error(f"Problem parsing entry {primary_AC}: failed to parse AC/OX/GN, debug me!")
            raise Exception('failed to parse uniprot entry')


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
