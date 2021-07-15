import argparse
import os
import subprocess
import re

"""
Script made by Javier Martín de Benito while businesss practicing at Instituto de
Tecnológico Agrario de Castilla y León (Valladolid, Castile and Leon).

"""

AUTHOR = "Javier Martín"
URL = "https://github.com/javiermdb99"

COLUMNS = ["FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
           "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", 
           "PRODUCT", "RESISTANCE"]
BLAST_FIELDS = ["qseqid", "qstart", "qend", "qlen", "ssequid", "sstart", "send",
                "slen", "sstrand", "evalue", "length", "pident", "gaps", "gapopen",
                "stitle"]
REQUIRE = ["blastn", "blastx", "makeblastdb", "blastdbcmd", "any2fasta", "gzip", "unzip"]

parser = argparse.ArgumentParser(description="A program")


def parseArguments():
    """
    Provides different options to change the program's behaviour
    """

    parser.add_argument('file', action='store', help=".fasta File.")
    parser.add_argument('-d', '--debug', action='store_true', help="Verbose debug output.")
    parser.add_argument('-q', '--quiet', action='store_true', help="No error output.")
    parser.add_argument('-v', '--version', action='store_true', help="Print version.")
    # parser.add_argument('--check', action='store_true', help="Check dependencies are installed.")
    parser.add_argument('-t', '--threads', action='store', help="Use this many BLAST+ threads.",
                        default=1, type=int)
    parser.add_argument('--fofn', action='store', help="Run files listed on this file.")
    parser.add_argument('-l', '--list', action='store_true', help="List included databases.")
    # parser.add_argument('--datadir', action='store_true', help="Verbose debug output") # TODO
    # parser.add_argument('--db', action='store_true', help="Verbose debug output")
    parser.add_argument('--noheader', action='store_true', help="Suppress column headers.")
    parser.add_argument('--csv', action='store_true', help="Output csv instead of tsv.")
    # parser.add_argument('--nopath', action='store_true', help="Strip filename paths from FILE column.")

    return vars(parser.parse_args())



def readFile(file):
    try:
        f = open(file, 'r')
    except IOError:
        print(file, " does not exist, or is unreadable.")
    # x=0
    # for _ in f:
    #     x = x+1
    
    print("Processing ", file)
    # print(x)

def blastDatabaseInfo():
    out = subprocess.Popen(['blastdbcmd', '-info', '-db', 
    '/home/javier/anaconda3/db/resfinder/sequences'], stdout=subprocess.PIPE)
    stdout, stderr = out.communicate()
    print(stdout)

args = parseArguments()
# print(args['fofn'])
# if args['threads'] < 1:
#     raise NameError('Threads must be greater or equal to 1.')
# # TODO CHEQUEOS
file = args['file']
readFile(file)
blastDatabaseInfo()
#os.system("echo Hola")
