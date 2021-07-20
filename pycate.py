import argparse
import sys
import subprocess
import re
from prettytable import PrettyTable
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
REQUIRE = ["blastn", "blastx", "makeblastdb",
           "blastdbcmd", "any2fasta", "gzip", "unzip"]


def process_debug(message):
    if debug:
        print(message)


def parse_arguments():
    """
    Provides different options to change the program's behaviour.
    """

    parser.add_argument('file', action='store', help=".fasta File.")
    parser.add_argument('-d', '--debug', action='store_true',
                        help="Verbose debug output.")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="No error output.")
    parser.add_argument('-v', '--version',
                        action='store_true', help="Print version.")
    # parser.add_argument('--check', action='store_true', help="Check dependencies are installed.")
    parser.add_argument('-t', '--threads', action='store',
                        help="Use this many BLAST+ threads.",
                        default=1, type=int)
    parser.add_argument('--fofn', action='store',
                        help="Run files listed on this file.")
    parser.add_argument('-l', '--list', action='store_true',
                        help="List included databases.")
    # parser.add_argument('--datadir', action='store_true', help="Verbose debug output") # TODO
    parser.add_argument('--db', action='store',
                        help="Database to be used.", default="resfinder")
    parser.add_argument('--noheader', action='store_true',
                        help="Suppress column headers.")
    parser.add_argument('--csv', action='store_true',
                        help="Output csv instead of tsv.")
    # parser.add_argument('--nopath', action='store_true', help="Strip filename paths from FILE column.")
    parser.add_argument('--minid', action='store', type=int, default=80,
                        help="Minimum DNA identity.")
    parser.add_argument('--mincov', action='store', type=int, default=80,
                        help="Minimum DNA coverage.")
    return vars(parser.parse_args())


def check_arguments(args):
    print(type(args['minid']))
    minid = args['minid']
    mincov = args['mincov']
    if not (minid <= 100 and minid > 0):
        raise ValueError('minid must be between 0 and 100 (100 included).')
    if not (mincov <= 100 and mincov >= 0):
        raise ValueError('mincov must be between 0 and 100 (both included).')


def blast_database_info(db):
    """
    This function parses the info output of the database that is going to be
    used.

    Arguments:
        -db: database used.

    Returns:
        -seq: number of sequences found in database.
        -total_bases: number of bases found in database.
        -date: date in which database was updated.
    """
    process_debug("Processing database information.")

    # Execute info command in BLAST and parse output
    out = subprocess.run(['blastdbcmd', '-info', '-db',
                          '/home/javier/anaconda3/db/' + db + '/sequences'], 
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    info = ' '.join(re.split("[ \n\t]", out.stdout.decode('utf-8')))

    try:
        seq = re.findall("[\d,]* sequences", info)[0]
        seq = int(seq.replace(",", "").replace(" sequences", ''))

        total_bases = re.findall("[\d,]* total bases", info)[0]
        total_bases = int(total_bases.replace(
            ",", '').replace(" total bases", ""))

        date = re.findall("Date: \w*\s+\d+,\s+\d+", info)[0]
        date = date.replace("Date: ", '')

        type = re.findall("total\s+\S+", info)[0]
        type = type.replace("total ", "")
        type = "prot" if type == "residues" else "nucl"

        # TODO
    except IndexError:
        print(db, "database does not exist.")
        sys.exit(1)

    return seq, total_bases, date, type


def process_file(file, type):
    process_debug("Reading file.")

    try:
        f = open(file, 'r')
    except IOError:
        print(file, " does not exist, or is unreadable.")

    print("Processing ", file)
    process_debug("Processing file. BLAST+ query")
    if type == 'nucl':
        blast_query = "blastn -task blastn -dust no -perc_identity " + \
            str(minid)
    else:
        blast_query = "blastx -task blastx-fast -seg no"

    query = ("(any2fasta -q -u " + file + " | " +
             blast_query + "-db \Q$db_path\E -outfmt '$format' -num_threads $threads" +
             " -evalue 1E-20 -culling_limit $CULL" +
             " -max_target_seqs 10000")

    # TODO
    #   my $blastcmd = $dbinfo->{DBTYPE} eq 'nucl'
    #            ? "blastn -task blastn -dust no -perc_identity $minid"
    #            : "blastx -task blastx-fast -seg no"
    #            ;
    # any2fasta -q -u MS7593.fasta | blastx -task blastx-fast -seg no -db db/resfinder/sequences -num_threads 1 -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
    # (any2fasta -q -u MS7593.fasta | blastn -task blastn -dust no -perc_identity 80  -db db/ncbi/sequences -num_threads 1 -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000 -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle')


parser = argparse.ArgumentParser(description="A program")
args = parse_arguments()
check_arguments(args)
file = args['file']
db = args['db']
debug = args['debug']
minid = args['minid']
mincov = args['mincov']


sequences, total_bases, date_db, type = blast_database_info(db)
print(f"\nDatabase \"{db}\" information:")
print(
    f"Number of sequences: {sequences}\t Number of bases: {total_bases}\t Date: {date_db}")
process_file(file, type)

t = PrettyTable(COLUMNS)

print(t)
