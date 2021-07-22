#!/usr/bin/env perl

use warnings;
use strict;
use File::Spec;
use Data::Dumper;
use FindBin;
use File::Spec;
use File::Slurp;
use List::MoreUtils qw(zip);
use Cwd qw(abs_path);
#use lib "$FindBin::RealBin/../perl5";

#..............................................................................
# Globals needed before --help etc

my $VERSION = "0.8";
my $EXE = just_filename( $FindBin::RealScript );
my $AUTHOR = 'Torsten Seemann';
my $EMAIL = 'torsten.seemann@gmail.com';
my $URL = 'https://github.com/tseemann/abricate';

#..............................................................................
# Command line options

my(@Options, $debug, $quiet, $setupdb, $list, $summary, $check,
             $datadir, $db, $minid, $mincov,
             $noheader, $csv, $nopath);
setOptions();

#..............................................................................
# Globals

my $OUTSEP = "\t";
my $ABSENT = '.';
my $FIELD = '%COVERAGE';
my $FIELDSEP = ';';
my $IDSEP = '~~~';

my $CULL = 1;   # BLAST -culling_limit 

my @BLAST_FIELDS = qw(
  qseqid qstart qend qlen
  sseqid sstart send slen sstrand
  evalue length pident gaps gapopen
  stitle
);

my @COLUMNS = qw(
 FILE SEQUENCE START END GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT
);
$COLUMNS[0] = '#'.$COLUMNS[0];

my @REQUIRE = qw(blastn makeblastdb blastdbcmd seqret gzip unzip);

#..............................................................................
# Option parsing

($minid > 0 and $minid <= 100) or err("--minid must be between 1 and 100");
($mincov >= 0 and $mincov <= 100) or err("--mincov must be between 0 and 100");
$OUTSEP = ',' if $csv;  # default is tab

if ($summary) {
  summary_table( @ARGV );
  exit;
}

if ($check) {
  $debug=1;
  msg("Checking dependencies are installed:");
  require_exe( @REQUIRE );
  msg("OK.");
  exit;
}

# check if dependencies are installed
require_exe( @REQUIRE );

# do these now that we know the deps are ok
if ($list or $setupdb) {
  list_databases($datadir, $setupdb);
  exit;
}

# check if blastn > 2.2.30 to support 'gaps' custom field
my($version) = qx(blastn -version 2>&1);
$version =~ m/2.(\d+.\d+)\+/ or err("Could not parse blastn version from '$version'");
$1 >= 2.30 or err("You need to install BLAST+ 2.2.30 or higher");

# check database exists
$datadir = abs_path($datadir);
my $db_path = "$datadir/$db/sequences";
my $dbinfo = blast_database_info($db_path);
$dbinfo or err("BLASTN database not found: $db_path");
msg("Using database $db: ", $dbinfo->{SEQUENCES}, "sequences - ", $dbinfo->{DATE});

# output format
my $format = "6 @BLAST_FIELDS";
print line(@COLUMNS) unless $noheader;

FILE:
for my $file (@ARGV) {
  my %seen;
  my @hit;
  msg("Processing: $file");
  -r $file or err("'$file' does not exist, or is unreadable");

  my $cmd = "(gzip -d -c -f \Q$file\E | seqret -auto -filter -osformat2 fasta |"
          . " blastn -db \Q$db_path\E -outfmt '$format'"
          . " -task blastn -evalue 1E-20 -dust no -max_target_seqs 10000 -perc_identity $minid"
          . " -culling_limit $CULL) 2>&1"
          ;
  msg("Running: $cmd") if $debug;
  
  open BLAST, "$cmd |";
  while (<BLAST>) { 
    chomp;
    my @x = split m/\t/;

    @x == @BLAST_FIELDS or err("WARNING: can not find sequence data in '$file'");

    my %hit = (map { $BLAST_FIELDS[$_] => $x[$_] } 0 .. @BLAST_FIELDS-1); 
    ($hit{sstart}, $hit{send}) = ($hit{send}, $hit{sstart}) if $hit{sstrand} eq 'minus';

    next if $seen{ join('~', @hit{qw(qseqid qstart qend)}) }++;
    
    my $pccov = 100 * ($hit{length}-$hit{gaps}) / $hit{slen};
    next unless $pccov >= $mincov;

#    $hit{sseqid} =~ s/_.*$//;
    my($database, $gene, $acc) = split m/$IDSEP/, $hit{sseqid};
    # if it wasn't in the format we expected, try and be sensible
    if (!defined $gene and !defined $acc) {
      $acc = '';
      $gene = $hit{sseqid};
      $database = $db;  # name specified with --db
    }
    msg( Dumper(\%hit) ) if $debug;

    my $product = $hit{'stitle'} || "n/a";
    $product =~ s/[,\t]//g;  # remove output separators

    my $minimap = minimap( @hit{qw(sstart send slen gapopen)} );
#    $minimap .= '!-'.$hit{gaps} if $hit{gaps} > 0;
    push @hit, [ 
      $nopath ? just_filename($file) : $file ,
      $hit{qseqid}, 
      $hit{qstart}, 
      $hit{qend},
      $gene,  #$hit{sseqid},
      $hit{sstart}.'-'.$hit{send}.'/'.$hit{slen}, 
      $minimap, 
      $hit{gapopen}.'/'.$hit{gaps},
      sprintf("%.2f", $pccov),
      sprintf("%.2f", $hit{pident}),
      $database,
      $acc,
      $product,
    ];
  }
  close BLAST;

  msg("Found", scalar(@hit), "genes in $file");

  # Sort hits
  #0    1   2     3   4    5
  #FILE CHR START END GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION
#  @hit = sort { $a->[4] cmp $b->[4] || $b->[8] <=> $a->[8] } @hit;
  @hit = sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @hit;

  print line(@$_) for @hit;  
}

#----------------------------------------------------------------------

sub minimap {
  my($x, $y, $L, $broken, $scale, $on, $off) = @_;
  my $WIDTH = 15 - ($broken ? 1 : 0);
  $broken ||= 0;
  $scale ||= ($L/$WIDTH);
  $on  ||= '=';
  $off ||= '.';
  $x = int( $x / $scale );
  $y = int( $y / $scale );
  $L = int( $L / $scale );
  my $map='';
  for my $i (0 .. $WIDTH-1) {
#    msg("$i $x $y $L $scale");
    $map .= ($i >= $x and $i <= $y) ? $on : $off;
    $map .= '/' if $broken and $i==int($WIDTH/2);
  }
  return $map;
}

#----------------------------------------------------------------------

sub summary_table {
  my(@fname) = @_;
  
  my %gene;  # genes observed
#  my %data = ( map { $_ => {} } @fname); 
  my %data;
  my @hdr;
  my %id_of;

  for my $fname (@fname) {
    err("Already summarized file '$fname'") if exists $data{$fname};
    # initialise to empty so it appears in report even if 0 genes found
    $data{$fname} = {};
    msg("Parsing: $fname") if $debug;
    open my $fh, '<', $fname or msg("WARNING: Can not open '$fname' to summarize.");
    while (<$fh>) {
      chomp;
      my @col = split m/$OUTSEP/;
      print STDERR "### [$.] @col\n" if $debug;
      @hdr = @col if ! @hdr;  # first row we see will be header
      next if $col[0] =~ m/^#/;
      $col[0] = just_filename($col[0]) if $nopath;
      $id_of{$fname} ||= $col[0];  # keep alternate "original filename"
      my $row = { zip @hdr, @col };
      print STDERR Dumper($row) if $debug;
      $gene{ $row->{'GENE'} }++;
      push @{ $data{ $fname }{ $row->{'GENE'} } } , $row;
    }
  }

  my @gene = sort { $a cmp $b } keys %gene;
  print line('#FILE', 'NUM_FOUND', @gene);

  for my $fname (sort { $a cmp $b } keys %data) {
    my @present = map { exists $data{$fname}{$_} 
                        ? join( $FIELDSEP, map { $_->{$FIELD} } @{$data{$fname}{$_}} ) 
                        : $ABSENT 
                      } @gene;
    print line( $fname || $fname, scalar( keys %{$data{$fname}} ), @present );
    # print line( $id_of{$fname} || $fname, scalar( keys %{$data{$fname}} ), @present );
  }
}

#----------------------------------------------------------------------

sub blast_database_info {
  my($prefix) = @_;
  my(@res) = qx(blastdbcmd -info -db \Q$prefix\E 2> /dev/null);
  chomp @res;
  my $res = join(' ', @res);
  my $info = { PREFIX => $prefix };

  $res =~ m/\b([\d,]+)\s+sequences;/ or return;
  $info->{SEQUENCES} = $1;
  $info->{SEQUENCES} =~ s/,//g;

  $res =~ m/\bDatabase: (\S+)/;
  $info->{TITLE} = $1;
  
  $res =~ m/\bDate: (\w+)\s+(\d+),\s+(\d{4})\s/;
  $info->{DATE} = "$3-$1-$2";  # YYYY-MM-DD
  
  return $info;
}

#----------------------------------------------------------------------

sub list_databases {
  my($from, $setup_too) = @_;
  my @dir = grep { -d "$from/$_" } read_dir($from);
  my $count=0;
  my @list;
  
  for my $name (sort @dir) {
    my $dbname = "$from/$name/sequences";
    if (-r $dbname) {
      if (not -r "$dbname.nin") {
        if ($setup_too) {
          msg("Database $name has not been indexed; formatting now.");
          system("makeblastdb -in '$dbname' -title '$name' -dbtype nucl -parse_seqids -hash_index -logfile $from/$name/makeblastdb.log");
        }
        else {
          err("Database $name is not indexed, please try: abricate --setupdb");
        }          
      }
      my $info = blast_database_info( $dbname );
      # save for end
      push @list, [ $name, $info->{SEQUENCES}, $info->{DATE} ];
    }
  }
  if (@list) {
    # print to STDOUT
    print line("DATABASE", "SEQUENCES", "DATE");
    print map { line(@$_) } @list;
  }
  else {
    msg("No databases found in $from");
  }
}

#----------------------------------------------------------------------

sub line {
  return join($OUTSEP, @_)."\n";
}

#----------------------------------------------------------------------

sub version {
  print "$EXE $VERSION\n";
  exit;
}

#----------------------------------------------------------------------

sub msg {
  print STDERR "@_\n" unless $quiet;
}

#----------------------------------------------------------------------

sub err {
  print STDERR "ERROR: @_\n";
  exit(1);
}

#----------------------------------------------------------------------

sub just_filename {
  my($path) = @_;
  my(undef,undef,$fname) = File::Spec->splitpath( $path );
  return $fname;
}

#----------------------------------------------------------------------

sub require_exe {
  my(@arg) = @_;
  for my $exe (@arg) {
    my $where = '';
    for my $dir ( File::Spec->path ) {
      if (-x "$dir/$exe") {
        $where = "$dir/$exe";
        last;
      }
    }
    if ($where) {
      msg("Found '$exe' => $where") if $debug;
    }
    else {
      err("Could not find '$exe'. Please install it and ensure it is in the PATH.");
    }
  }
  return;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Verbose debug output"},
    {OPT=>"quiet!",  VAR=>\$quiet, DEFAULT=>0, DESC=>"Quiet mode, no stderr output"},
    {OPT=>"version!",  VAR=>\&version, DESC=>"Print version and exit"},
    {OPT=>"setupdb!",  VAR=>\$setupdb, DEFAULT=>0, DESC=>"Format all the BLAST databases"},
    {OPT=>"list!",  VAR=>\$list, DEFAULT=>0, DESC=>"List included databases"},
    {OPT=>"check!",  VAR=>\$check, DEFAULT=>0, DESC=>"Check dependencies are installed"},
    {OPT=>"summary!",  VAR=>\$summary, DEFAULT=>0, DESC=>"Summarize multiple reports into a table"},
    {OPT=>"datadir=s",  VAR=>\$datadir, DEFAULT=>"$FindBin::RealBin/../db", DESC=>"Location of database folders"},
    {OPT=>"db=s",  VAR=>\$db, DEFAULT=>"resfinder", DESC=>"Database to use"},
    {OPT=>"noheader!",  VAR=>\$noheader, DEFAULT=>0, DESC=>"Suppress column header row"},
    {OPT=>"csv!",  VAR=>\$csv, DEFAULT=>0, DESC=>"Output CSV instead of TSV"},
    {OPT=>"minid=f",  VAR=>\$minid, DEFAULT=>75, DESC=>"Minimum DNA %identity"},
    {OPT=>"mincov=f",  VAR=>\$mincov, DEFAULT=>0, DESC=>"Minimum DNA %coverage"},
    {OPT=>"nopath!",  VAR=>\$nopath, DEFAULT=>0, DESC=>"Strip filename paths from FILE column"},
  );

  @ARGV or usage(1);

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage(1);

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  my($exitcode) = @_;
  $exitcode = 0 if $exitcode eq 'help'; # what gets passed by getopt func ref
  $exitcode ||= 0;
  select STDERR if $exitcode; # write to STDERR if exitcode is error

  print "Synopsis:\n\tFind and collate amplicons in assembled contigs\n";
  print "Author:\n\t$AUTHOR <$EMAIL>\n";
  print "Usage:\n";
  print "\t% $EXE --list\n";
  print "\t% $EXE [options] <contigs.{fasta,gbk,embl}[.gz]> > out.tab\n";
  print "\t% $EXE --summary <out1.tab> <out2.tab> <out3.tab> ... > summary.tab\n";
  print "Options:\n";
  foreach (@Options) {
    my $opt = $_->{OPT};
    $opt =~ s/!$//;
    $opt =~ s/=s$/ [X]/;
    $opt =~ s/=i$/ [N]/;
    $opt =~ s/=f$/ [n.n]/;
    printf "\t--%-13s %s%s.\n",$opt,$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  print "Documentation:\n\t$URL\n";
  exit($exitcode);
}
 
#----------------------------------------------------------------------