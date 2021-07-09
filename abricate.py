import argparse

"""
Script made by Javier Martín de Benito while businesss practicing at Instituto de
Tecnológico Agrario de Castilla y León (Valladolid, Castile and Leon).

"""

__author__ = "Javier Martín"
__url__ = "https://github.com/javiermdb99"

columns = "FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE".split()

def parseArguments():

    """
    Provides different options to change the program's behaviour
    """
    parser = argparse.ArgumentParser(description="A program")

    parser.add_argument('-d', '--debug', action='store_true', help="Verbose debug output.")
    parser.add_argument('-q', '--quiet', action='store_true', help="No error output.")
    parser.add_argument('-v', '--version', action='store_true', help="Print version and exit.")
    # parser.add_argument('--check', action='store_true', help="Check dependencies are installed.")
    parser.add_argument('-t', '--threads', action='store', help="Use this many BLAST+ threads.",
                        default=1, type=int)
    # parser.add_argument('--fofn', action='store_true', help="Run files listed on this file.") # TODO leer fichero
    parser.add_argument('-l', '--list', action='store_true', help="List included databases.")
    # parser.add_argument('--datadir', action='store_true', help="Verbose debug output") # TODO
    # parser.add_argument('--db', action='store_true', help="Verbose debug output")
    parser.add_argument('--noheader', action='store_true', help="Suppress column headers.")
    parser.add_argument('--csv', action='store_true', help="Output csv instead of tsv.")
    # parser.add_argument('--nopath', action='store_true', help="Strip filename paths from FILE column.")

    parser.parse_args()
    # TODO filtering

    # TODO CHEQUEOS

# TODO:
#   Realizar la función de lectura del fasta
#   Poder añadirla como un argumento necesario en el parseArguments.
def readFile(file):
    f = open(file, 'r')
    i = 0
    for _ in f:
        i+=1
    print(i)

parseArguments()