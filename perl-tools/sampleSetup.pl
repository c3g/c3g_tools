#!/usr/bin/perl

use Cwd;
use File::Basename;
use File::Path qw(mkpath);
use Getopt::Long;
#use Text::CSV;
use Text::CSV::Encoded;

my $version = "1.0";

use strict;

&main();

sub getUsage {
  my $usage = <<END;
Usage: perl $0 --nanuqAuthFile \$HOME/.nanuqAuth.txt --usesheet project.nanuq.csv --tech HiSeq
  --nanuqAuthFile  <FILE>          Path to Nanuq authentication file
  --projectId      <INT>           Nanuq project ID from which to get the sample sheet (can't be used with --usesheet)
  --usesheet       <FILE>          Use the specified sample sheet instead of fecthing it from Nanuq (can't be used with --projectId)
  --format         [bam|fastq]     Raw reads format (default: bam)
  --links                          Create raw_reads directory and symlinks (default)
  --nolinks                        Do not create raw_reads directory or symlinks
  --tech           [HiSeq|MiSeq]   Sequencing technology
  --help                           Show this help

  The 'nanuqAuthFile' contains your Nanuq username and password.
  To create it, type the command:
  echo -n "user=<USERNAME>&password=<PASSWD>" > \$HOME/.nanuqAuth.txt ; chmod u+r,go-rwx \$HOME/.nanuqAuth.txt
  The '-n' is important because there cannot be a Carriage Return, EOL at the end of the line.
END

  return $usage;
}

sub main {
  my $techName;
  my $projectId;
  my $sampleSheet;
  my $format = "bam";  # Raw reads are in bam format by default, otherwise in fastq format if specified
  my $nanuqAuthFile;
  my $links = 1;  # Create symlinks by default
  my $help;
  my $result = GetOptions(
    "tech=s"          => \$techName,
    "projectId=i"     => \$projectId,
    "usesheet=s"      => \$sampleSheet,
    "format=s"        => \$format,
    "nanuqAuthFile=s" => \$nanuqAuthFile,
    "links!"          => \$links,
    "help!"           => \$help,
  );

  if ($help) {
    die getUsage();
  }

  my $errMsg = "";
  if (!defined($nanuqAuthFile) or !-e $nanuqAuthFile) {
    $errMsg .= "Error: missing nanuqAuthFile!\n";
  }
  if (defined($projectId) and defined($sampleSheet)) {
    $errMsg .= "Error: --projectId and --useSheet options cannot be used together!\n";
  }
  if ((!defined($projectId) or length($projectId) == 0) and (!defined($sampleSheet) or length($sampleSheet) == 0)) {
    $errMsg .= "Error: missing --projectId or --useSheet option!\n";
  }
  if ($format ne "bam" and $format ne "fastq") {
    $errMsg .= "Error: invalid --format value (should be 'bam' or 'fastq')!\n";
  }
  if (!defined($techName) or length($techName) == 0 or ($techName ne "HiSeq" and $techName ne "MiSeq")) {
    $errMsg .= "Error: missing or invalid --tech value (should be 'HiSeq' or 'MiSeq')!\n";
  }
  if (length($errMsg)) {
    die $errMsg . "\n" . getUsage();
  }

  # Default Nanuq project file name
  my $projectFile = 'project.nanuq.csv';

  if (defined($projectId)) {
    # Fecth sample sheet data from Nanuq
    getSampleSheet($projectFile, $techName, $projectId, $nanuqAuthFile);
  } else {
    # Use the specified sample sheet
    $projectFile = $sampleSheet;
  }

  if ($links) {
    my $rA_sampleInfos = parseSampleSheet($projectFile, $format, $techName);
    createLinks($rA_sampleInfos);
  }
}

sub getSampleSheet {
  my $projectFile = shift;
  my $techName = shift;
  my $projectId = shift;
  my $nanuqAuthFile = shift;

  my $command = 'wget --no-cookies --post-file ' . $nanuqAuthFile . ' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/' . $techName . '/project/' . $projectId . '/filename/' . $projectFile . "\n";
  print '#' . $command;
  system($command);
  if ($? == -1) {
    print "failed to execute: $!\n";
    exit(1);
  } elsif ($? & 127) {
    printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
    exit(1);
  } else {
    my $childValue = $? >> 8;
    if ($childValue != 0) {
      printf "child exited with value %d\n", $childValue;
      exit(1);
    }
  }
}

sub parseSampleSheet {
  my $fileName = shift;
  my $format = shift;
  my $techName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $readSetIdIdx=-1;
  my $filePrefixIdx=-1;

  my $csv = Text::CSV::Encoded->new ({ encoding => "iso-8859-1" });
  $csv->parse($line);
  my @headers = $csv->fields();
  for (my $idx = 0; $idx < @headers; $idx++) {

    $headers[$idx] =~ s/"//g;
    if ($headers[$idx] eq "Name") {
      $nameIdx = $idx;
    } elsif ($headers[$idx] eq "Library Barcode") {
      $libraryBarcodeIdx = $idx;
    } elsif ($headers[$idx] eq "Run") {
      $runIdIdx = $idx;
    } elsif ($headers[$idx] eq "Region") {
      $laneIdx = $idx;
    } elsif ($headers[$idx] eq "Run Type") {
      $runTypeIdx = $idx;
    } elsif ($headers[$idx] eq "Status") {
      $statusIdx = $idx;
    } elsif ($headers[$idx] eq "Read Set Id") {
      $readSetIdIdx = $idx;
    } elsif ($headers[$idx] eq "Filename Prefix") {
      $filePrefixIdx = $idx;
    }
  }

  my $sampleSheetErrors = "";
  if ($nameIdx == -1) {
    $sampleSheetErrors .= "Missing Sample Name\n";
  }
  if ($libraryBarcodeIdx == -1) {
    $sampleSheetErrors .= "Missing Library Barcode\n";
  }
  if ($runIdIdx == -1) {
    $sampleSheetErrors .= "Missing Run ID\n";
  }
  if ($laneIdx == -1) {
    $sampleSheetErrors .= "Missing Lane\n";
  }
  if ($runTypeIdx == -1) {
    $sampleSheetErrors .= "Missing Run Type\n";
  }
  if ($statusIdx == -1) {
    $sampleSheetErrors .= "Missing Status\n";
  }
  if ($readSetIdIdx == -1) {
    $sampleSheetErrors .= "Missing Read Set Id\n";
  }
  if ($filePrefixIdx == -1) {
    $sampleSheetErrors .= "Missing Filename Prefix\n";
  }
  if (length($sampleSheetErrors) > 0) {
    die $sampleSheetErrors;
  }

  while ($line = <SAMPLE_SHEET>) {
    $csv->parse($line);
    my @values = $csv->fields();
    if ($values[$statusIdx] =~ /invalid/) {
      warn "[Warning] Sample Name $values[$nameIdx], Run ID $values[$runIdIdx], Lane $values[$laneIdx] data is invalid!\n";
    } else {
      my %sampleInfo;
      $sampleInfo{'name'} = $values[$nameIdx];
      $sampleInfo{'libraryBarcode'} = $values[$libraryBarcodeIdx];
      $sampleInfo{'runId'} = $values[$runIdIdx];
      $sampleInfo{'lane'} = $values[$laneIdx];
      $sampleInfo{'runType'} = $values[$runTypeIdx];
      $sampleInfo{'readSetId'} = $values[$readSetIdIdx];
      $sampleInfo{'filePrefix'} = $values[$filePrefixIdx];

      my $rootDir;
      if ($techName eq 'HiSeq') {
        $rootDir = "/lb/robot/hiSeqSequencer/hiSeqRuns/";
      } elsif ($techName eq 'MiSeq') {
        $rootDir = "/lb/robot/miSeqSequencer/miSeqRuns/";
      } else {
        die "Unknown prefix technology type: " . $sampleInfo{'filePrefix'} . "\n";
      }

      # Find yearly directories
      opendir(ROOT_DIR, $rootDir) or die "Couldn't open directory " . $rootDir . "\n";
      my @roots = grep { /^2\d\d\d/ } readdir(ROOT_DIR);
      closedir(ROOT_DIR);
      for (my $i = 0; $i < @roots; $i++) {
        $roots[$i] = $rootDir . '/' . $roots[$i];
      }
      push(@roots, $rootDir);

      my @rootFiles;
      my $runPath;
      for my $rootDir (@roots) {
        my @tmpPaths;
        opendir(ROOT_DIR, $rootDir) or die "Couldn't open directory " . $rootDir . "\n";
        if ($techName eq 'HiSeq') {
          @tmpPaths = grep { /.*[0-9]+_[^_]+_[^_]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR);
        } else {
          @tmpPaths = grep { /.*[0-9]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR);
        }

        if (@tmpPaths > 0) {
          push(@rootFiles, @tmpPaths);
          $runPath = $rootDir . '/' . $rootFiles[0];
        }
      }

      if (@rootFiles == 0) {
        die "Run not found: " . $sampleInfo{'runId'} . "\n";
      } elsif (@rootFiles > 1) {
        die "Many runs found: " . $sampleInfo{'runId'} . "\n";
      }

      my $rawReadDir = `echo $runPath/se*`;
      chomp($rawReadDir);

      my $rawReadFile1;
      my $rawReadFile2;

      if ($rawReadDir =~ /\*/) {
        $rawReadDir = `echo $runPath/Data/In*/B*/G*`;
        chomp($rawReadDir);
        if ($rawReadDir =~ /\*/) {
          die "Couldn't find fastq directory: $rawReadDir\n";
        }

        $sampleInfo{'qualOffset'} = "64";

        $rawReadFile1 = $rawReadDir . '/s_' . $sampleInfo{'lane'} . '_1_*' . $sampleInfo{'name'} . '*.txt.gz';
        $rawReadFile2 = $rawReadDir . '/s_' . $sampleInfo{'lane'} . '_2_*' . $sampleInfo{'name'} . '*.txt.gz';
        $sampleInfo{'filename2'} = `echo $rawReadFile2`;
      } else {
        $sampleInfo{'qualOffset'} = "33";

        if ($format eq "fastq") {
          $rawReadFile1 = $rawReadDir . '/' . $sampleInfo{'filePrefix'} . '_R1.fastq.gz';

          my $runType = $values[$runTypeIdx];
          if ($runType eq "PAIRED_END") {
            $rawReadFile2 = $rawReadDir . '/' . $sampleInfo{'filePrefix'} . '_R2.fastq.gz';
            $sampleInfo{'filename2'} = `echo $rawReadFile2`;
          }
        } else {    # Bam format by default
          $rawReadFile1 = $rawReadDir . '/' . $sampleInfo{'filePrefix'} . '.bam';
        }
      }

      $sampleInfo{'filename1'} = `echo $rawReadFile1`;
      chomp($sampleInfo{'filename1'});
      if ($values[$runTypeIdx] eq "PAIRED_END") {
        chomp($sampleInfo{'filename2'});
      }
      push(@retVal, \%sampleInfo);
    }
  }

  return \@retVal;
}

sub createLinks {
  my $rA_sampleInfos = shift;

  for my $rH_sample (@$rA_sampleInfos) {
    my $rawReadDir = 'raw_reads/' . $rH_sample->{'name'} . "/run" . $rH_sample->{'runId'} . "_" . $rH_sample->{'lane'};
    mkpath($rawReadDir);

    my $rawReadPrefix = $rawReadDir . '/' . $rH_sample->{'name'} . '.' . $rH_sample->{'libraryBarcode'} . '.' . $rH_sample->{'qualOffset'} . ".";

    my @symlinks;

    # Retrieve target file extension and list all links to create
    if ($rH_sample->{'filename1'} =~ /\.bam$/) {    # BAM format
      push(@symlinks, [$rH_sample->{'filename1'}, $rawReadPrefix . "bam"]);
    } else {    # FASTQ format
      my $runType = $rH_sample->{'runType'};
      if ($runType eq "SINGLE_END") {
        push(@symlinks, [$rH_sample->{'filename1'}, $rawReadPrefix . "single.fastq.gz"]);
      } elsif ($runType eq "PAIRED_END") {
        push(@symlinks, [$rH_sample->{'filename1'}, $rawReadPrefix . "pair1.fastq.gz"]);
        push(@symlinks, [$rH_sample->{'filename2'}, $rawReadPrefix . "pair2.fastq.gz"]);
      } else {
        die "Error: unknown run type: $runType!";
      }
    }

    # Create all symbolic links
    for my $symlink (@symlinks) {
      if (-l @$symlink[1]) {
        warn "[Warning] Symbolic link @$symlink[1] already exists! Skipping.\n";
      } elsif (-f @$symlink[0] and symlink(@$symlink[0], @$symlink[1])) {
        print "Created symbolic link @$symlink[1] successfully.\n";
      } else {
        die "[Error] Can't create symbolic link @$symlink[1] to target @$symlink[0]!\n";
      }
    }
  }
}
