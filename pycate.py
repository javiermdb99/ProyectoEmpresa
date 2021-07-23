import argparse
import sys
import subprocess
import re
import os
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
BLAST_FIELDS = ["qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send",
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
    # parser.add_argument('-q', '--quiet', action='store_true',
    #                     help="No error output.")
    # parser.add_argument('-v', '--version',
    #                     action='store_true', help="Print version.")
    # parser.add_argument('--check', action='store_true', help="Check dependencies are installed.")
    parser.add_argument('--threads', action='store',
                        help="Use this many BLAST+ threads.",
                        default=1, type=int)
    # parser.add_argument('--setupdb', action='store',
    #                     help="Formal all the BLAST databases.")
    # parser.add_argument('-l', '--list', action='store_true',
    #                     help="List included databases.")
    parser.add_argument('--datadir', action="store", default=os.getcwd() + "/db",
                        help="Directory where data is stored.")
    parser.add_argument('--db', action='store',
                        help="Database to be used.", default="resfinder")
    # parser.add_argument('--noheader', action='store_true',
    #                     help="Suppress column headers.")
    # parser.add_argument('--csv', action='store_true',
    #                     help="Output csv instead of tsv.")
    # parser.add_argument('--nopath', action='store_true', help="Strip filename paths from FILE column.")
    parser.add_argument('--minid', action='store', type=int, default=80,
                        help="Minimum DNA identity.")
    parser.add_argument('--mincov', action='store', type=int, default=0,
                        help="Minimum DNA coverage.")
    # parser.add_argument('--summary', action='store_true',
    #                     help="Summarize reports.")
    return vars(parser.parse_args())


def check_arguments(args):
    minid = args['minid']
    mincov = args['mincov']
    threads = args['threads']
    if not (minid <= 100 and minid > 0):
        raise ValueError('minid must be between 0 and 100 (100 included).')
    if not (mincov <= 100 and mincov >= 0):
        raise ValueError('mincov must be between 0 and 100 (both included).')
    if not (threads >= 1):
        raise ValueError('threads number must be greater than 0.')


def blast_database_info(db, db_name) -> tuple:
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
                          db],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    info = ' '.join(re.split("[ \n\t]", out.stdout.decode('utf-8')))

    try:
        seq = re.search("([\d,]*) sequences", info).group(1)

        total_bases = re.search("([\d,]*) total bases", info).group(1)

        date = re.search("Date: (\w*\s+\d+,\s+\d+)", info).group(1)

        type = re.search("total\s+(\S+)", info).group(1)
        type = "prot" if type == "residues" else "nucl"

    except:
        print(db_name, "database does not exist.")
        sys.exit(1)

    return seq, total_bases, date, type


def process_file(file, t: PrettyTable, type, db, db_name, threads, minid, mincov):
    process_debug("Reading file.")
    format = "6 " + " ".join(BLAST_FIELDS)

    print("Processing ", file)
    process_debug("Processing file. BLAST+ query")
    if type == 'nucl':
        blast_query = "blastn -task blastn -dust no -perc_identity " + \
            str(minid)
    else:
        blast_query = "blastx -task blastx-fast -seg no"

    query = (f"(any2fasta -q -u {file} | " +
             f" {blast_query} -db {db} -outfmt \'{format}\' -num_threads {threads}" +
             f" -evalue 1E-20 -culling_limit 1" +  # $CULL de abricate es el 1 aquí
             f" -max_target_seqs 10000)")
    cmd = subprocess.run(query,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         shell=True)
    output = cmd.stdout.decode('utf-8')
    output = output.split("\n")
    output.pop()

    print("Found ", len(output), " genes in ", file)
    # try:
    for i in output:
        line = i.split()
        output_dict = dict(zip(BLAST_FIELDS, line))
        if output_dict['sstrand'] == 'minus':
            (output_dict['sstart'], output_dict['send']) = (
                output_dict['send'], output_dict['sstart'])
        gene, acc = re.search("~{3}(\S+)~{3}(\S+)", output_dict['sseqid']).groups()
        row_cov = (100 * (int(output_dict['length']) - int(output_dict['gaps'])) /
                    int(output_dict['slen']))
        if row_cov < mincov :
            continue
        row_output = [file,
                        output_dict['qseqid'],
                        output_dict['qstart'],
                        output_dict['qend'],
                        # '-' if output_dict['sstrand'] == 'minus' else '+',
                        # output_dict['sstrand'],
                        gene,
                        f"{output_dict['sstart']}-{output_dict['send']}" +
                        f"/{output_dict['slen']}",
                        # minimap,
                        f"{output_dict['gapopen']}/{output_dict['gaps']}",
                        "{:.2f}".format(row_cov),
                        "{:.2f}".format(float(output_dict['pident'])),
                        db_name,
                        acc
                        ]

        t.add_row(row_output)
    # except:
    #     print("Sequence not found in ", file)
    #     sys.exit(1)

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
wd = os.path.abspath(args['datadir'])
threads = args['threads']

db_name = db
db = f"{wd}/{db}/sequences"

process_debug("Using " + wd + " working directory.")
process_debug("Using " + db + " database.")


sequences, total_bases, date_db, type = blast_database_info(db, db_name)
print(f"\nDatabase \"{db_name}\" information:")
print(
    f"Number of sequences: {sequences}\t Number of bases: {total_bases}\t Date: {date_db}")

t = PrettyTable(COLUMNS)
t = PrettyTable(['FILE', 'SEQUENCE', 'START',
                'END', 'GENE', 'COVERAGE', 'GAPS', '%COVERAGE', '%IDENTITY',
                'DATABASE', 'ACCESSION'])

#     COLUMNS = ["FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
#    "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION",
#    "PRODUCT", "RESISTANCE"]
process_file(file, t, type, db, db_name, threads, minid, mincov)

print(t)
