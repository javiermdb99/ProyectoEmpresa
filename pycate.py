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

# COLUMNS = ["FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
#            "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION",
#            "PRODUCT", "RESISTANCE"]

COLUMNS = ['FILE', 'SEQUENCE', 'START',
           'END', 'STRAND', 'GENE', 'COVERAGE', 'COVERAGE_MAP', 'GAPS', '%COVERAGE', 
           '%IDENTITY', 'DATABASE', 'ACCESSION']

BLAST_FIELDS = ["qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send",
                "slen", "sstrand", "evalue", "length", "pident", "gaps", "gapopen",
                "stitle"]
REQUIRE = ["blastn", "blastx", "makeblastdb",
           "blastdbcmd", "any2fasta", "gzip", "unzip"]


def process_debug(message):
    """
    Prints debug messages if needed.
    """
    if debug:
        print(message)


def parse_arguments():
    """
    Provides different options to change the program's behaviour.
    """

    # NOTE: OBSERVAR COMPORTAMIENTO DE NARGS
    parser.add_argument('file', action='store', nargs="?", help=".fasta File.")
    parser.add_argument('-d', '--debug', action='store_true',
                        help="Verbose debug output.")
    # parser.add_argument('-q', '--quiet', action='store_true',
    #                     help="No error output.")
    # parser.add_argument('--check', action='store_true', help="Check dependencies are installed.")
    parser.add_argument('--threads', action='store',
                        help="Use this many BLAST+ threads.",
                        default=1, type=int)
    parser.add_argument('--fofn', action='store_true',
                        help="Use file of file names.")
    parser.add_argument('--setupdb', action='store_true',
                        help="Formal all the BLAST databases.")
    parser.add_argument('-l', '--list', action='store_true',
                        help="List included databases.")
    parser.add_argument('--datadir', action="store", default=os.getcwd() + "/db",
                        help="Directory where data is stored.")
    parser.add_argument('--db', action='store',
                        help="Database to be used.", default="resfinder")
    parser.add_argument('--noheader', action='store_true',
                        help="Suppress column headers (this activates csv).")
    parser.add_argument('--csv', action='store_true',
                        help="Output csv instead of table.")
    parser.add_argument('--nopath', action='store_true',
                        help="Strip filename paths from FILE column.")
    parser.add_argument('--minid', action='store', type=int, default=80,
                        help="Minimum DNA identity.")
    parser.add_argument('--mincov', action='store', type=int, default=0,
                        help="Minimum DNA coverage.")
    # parser.add_argument('--summary', action='store_true',
    #                     help="Summarize reports.")
    return vars(parser.parse_args())


def check_arguments(args):
    """
    Checks all arguments are between allowed limits.

    Arguments:
        args: list of ParseArgument arguments
    """

    file = args['file']
    minid = args['minid']
    mincov = args['mincov']
    threads = args['threads']
    if file == None:
        raise ValueError('Please, specify a file name.')
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
        -db_name: name of database

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
        # Sequences found in database
        seq = re.search("([\d,]*) sequences", info).group(1)

        total_bases = re.search("([\d,]*) total bases", info).group(1)

        date = re.search("Date: (\w*\s+\d+,\s+\d+)", info).group(1)

        # nucleotides or proteines
        type = re.search("total\s+(\S+)", info).group(1)
        type = "prot" if type == "residues" else "nucl"

    except:
        print(db_name, "database does not exist.")
        sys.exit(1)

    return seq, total_bases, date, type

def gen_map(start, end, length, gaps=0):
    """
    Creates a simple representation of gene coverage.

    Arguments:
        -start: base number in which coverage starts.
        -end: base number in which coverage ends.
        -length: length of the gene.
        -gaps = if gene coverage is broken, this is greater than 0.
        Default value is 0.
    
    Returns:
        -gen_map
    """

    width = 15 - (1 if gaps else 0)
    scale = length // width
    on = '='
    off = '.'
    start = start // scale
    end = end // scale
    length = length // scale
    # print((gaps, width, start, end, width // 2, int(width/2)))
    gen_map = ''
    for i in range(1, width+1):
        current_gen = on if (i >= start and i <= end) else off
        gen_map = gen_map + current_gen
        if (gaps and i == (width // 2)):
            gen_map = gen_map + '/'
    return gen_map

def process_file(file, type, db, db_name, threads, minid, mincov, csv, no_path):
    """
    This function calls BLAST to process the file. Before doing this, it converts
    it into an allowed fasta file. Depending on which format it is being used, it
    will prompt it out as a table or a csv file.

    Arguments:
        -file: name of file
        -type: type of database being used (nucleotides or proteines)
        -db: route of database
        -db_name: name of database
        -threads: number of threads os BLAST query.
        -minid: minimum DNA identity.
        -mincov: minimum DNA coverage %.
        -csv: in case it is wanted to store the output in csv format.
    """

    file_name = os.path.basename(file) if no_path else file
    process_debug(f"Reading file {file}.")
    format = "6 " + " ".join(BLAST_FIELDS)

    if not csv:
        print("Processing ", file)

    # BLAST command
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

    if not csv:
        print("Found ", len(output), " genes in ", file)
    # try:
    for i in output:
        line = i.split()
        output_dict = dict(zip(BLAST_FIELDS, line))
        if output_dict['sstrand'] == 'minus':
            (output_dict['sstart'], output_dict['send']) = (
                output_dict['send'], output_dict['sstart'])
        gene, acc = re.search("~{3}(\S+)~{3}(\S+)",
                              output_dict['sseqid']).groups()
        row_cov = (100 * (int(output_dict['length']) - int(output_dict['gaps'])) /
                   int(output_dict['slen']))
        if row_cov < mincov:
            continue
        map = gen_map(int(output_dict["sstart"]), int(output_dict['send']),
                        int(output_dict["slen"]), int(output_dict["gapopen"]))
        row_output = [file_name,
                      output_dict['qseqid'],
                      output_dict['qstart'],
                      output_dict['qend'],
                      '-' if output_dict['sstrand'] == 'minus' else '+',
                      gene,
                      f"{output_dict['sstart']}-{output_dict['send']}" +
                      f"/{output_dict['slen']}",
                      map,
                      f"{output_dict['gapopen']}/{output_dict['gaps']}",
                      "{:.2f}".format(row_cov),
                      "{:.2f}".format(float(output_dict['pident'])),
                      db_name,
                      acc
                      ]
        if not csv:
            t.add_row(row_output)
        else:
            print(",".join(row_output))

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

def setup_database(path, name):
    # print(path)
    # print(name)
    try:
        open(path, 'r')
        execution = f"egrep '^>' {path}"
        cmd = subprocess.run(execution,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         shell=True)
        output = cmd.stdout.decode('utf-8')
        output = output.split('\n')
        output.pop()
        output = ''.join(output)
        print(output)
    except:
        print("Could not open ", name)

def list_databases(wd, t, setupdb):
    """
    This function searchs indexed databases in the given working directory.

    Arguments:
        - wd: working directory
        - t: PrettyTable which rows will be indexed databases information.
        - setupdb: if user wants to create a new database
    """

    subdirectories = [x[0] for x in os.walk(wd)]
    print(subdirectories)
    for subdirectory in subdirectories:
        name = re.search("/(\w+$)", subdirectory).group(1)
        db_name = f"{subdirectory}/sequences"
        try:
            open(db_name, 'r')
            # TODO: creación de la base de datos
            if setupdb:
                setup_database(db_name, name)

            # TODO: indexación de la base de datos con .nin
        except:
            continue
        seq, total_bases, date, type = blast_database_info(db_name, name)
        t.add_row([name, seq, date, type])
    return t


parser = argparse.ArgumentParser(description="A program")
args = parse_arguments()
wd = os.path.abspath(args['datadir'])
debug = args['debug']
list = args['list']
setupdb = args['setupdb']

# if setupdb:
#     setup_database("db/resfinder/sequences", "resfinder")
#     exit(0)

if list or setupdb:
    t = PrettyTable(["DATABASE", "SEQUENCE", "DATE", "TYPE"])
    t = list_databases(wd, t, setupdb)
    print(t)
else:
    check_arguments(args)
    file = args['file']
    db = args['db']
    minid = args['minid']
    mincov = args['mincov']
    threads = args['threads']
    fofn = args['fofn']
    no_header = args['noheader']
    csv = True if no_header else args['csv']
    no_path = args['nopath']

    db_name = db
    db = f"{wd}/{db}/sequences"

    process_debug("Using " + wd + " working directory.")
    process_debug("Using " + db + " database.")

    sequences, total_bases, date_db, type = blast_database_info(db, db_name)

    if not csv:
        print(f"\nDatabase \"{db_name}\" information:")
        print(
            f"Number of sequences: {sequences}\t Number of bases: {total_bases}\t Date: {date_db}")

    t = PrettyTable(COLUMNS)

    if csv and not no_header:
        print(",".join(COLUMNS))

    # in case fofn flag is used, file is where file names are stored
    if fofn:
        try:
            fofn_reader = open(file, 'r')
            files = fofn_reader.readlines()
            for file in files:
                process_file(file.strip(), type, db, db_name,
                             threads, minid, mincov, csv, no_path)
        except IOError:
            print("Can't open fofn ", file)
    else:
        process_file(file, type, db, db_name, threads,
                     minid, mincov, csv, no_path)

    if not csv:
        print(t)
