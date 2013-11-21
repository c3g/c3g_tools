#!/usr/bin/perl -w

=head1 NAME

I<Log>

=head1 SYNOPSIS

getLogReport <job_log_list_file>

=head1 DESCRIPTION

Parse internal log files and produce custom log reports

=head1 AUTHOR
B<Joel Fillon> - I<joel.fillon@mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

# Dependencies
#-----------------------
use POSIX qw(strftime);
use File::Basename;
use File::Spec;

# Check and assign the parameters
$ARGV[0] or die "Usage: perl " . basename($0) . " <job_log_list_file>\n";
my $jobLogList = $ARGV[0];

print getLogTextReport($jobLogList);

# SUB
#-----------------------

# Return an array of hashes of log values for each job's log file created by the cluster
sub readJobLogListFile {
  my $jobLogListPath = shift;
  # Return array of hash of log values
  my $rAoH_jobLogList = shift;

  # Read the global list of log files
  open(JOB_LOG_LIST_FILE, $jobLogListPath) or die "Cannot open $jobLogListPath\n";
  while(my $line = <JOB_LOG_LIST_FILE>) {
    chomp($line);

    # Retrieve each job's log file path
    my ($jobId, $jobName, $jobDependencies, $clusterJobLogPath) = split(/\t/, $line);

    if (defined $jobName and defined $clusterJobLogPath) {
      my %jobLog;

      $jobLog{'jobId'} = $jobId;
      $jobLog{'jobName'} = $jobName;
      $jobLog{'jobDependencies'} = $jobDependencies;

      # In old job log list version, cluster job log path was absolute. Now it is relative to job log list directory (useful if project directory has been moved or renamed).
      my $clusterJobLogFullPath = File::Spec->rel2abs(dirname($jobLogListPath)) . "/" . $clusterJobLogPath;
      unless (-e $clusterJobLogFullPath) {
        # Assume cluster job log path is absolute (old version)
        $clusterJobLogFullPath = $clusterJobLogPath;
      }
      $jobLog{'path'} = $clusterJobLogFullPath;

      # Read the job log file
      if (open(CLUSTER_JOB_LOG_FILE, $clusterJobLogFullPath)) {
        while(my $jobLine = <CLUSTER_JOB_LOG_FILE>) {
          # Job start date
          if ($jobLine =~ /^Begin PBS Prologue (.*) (\d+)$/) {
            $jobLog{'startDate'} = $1;
            $jobLog{'startSecondsSinceEpoch'} = $2;
          # Job number
          } elsif ($jobLine =~ /^Job ID:\s+(\S+)/) {
            $jobLog{'jobFullId'} = $1;
          # Job MUGQIC exit status
          } elsif ($jobLine =~ /MUGQICexitStatus:(\d+)/) {
            $jobLog{'MUGQICexitStatus'} = $1;
          # Job exit status (should be the same as MUGQIC exit status unless MUGQIC exit status is skipped)
          } elsif ($jobLine =~ /^Exit_status:\s+(\d+)/) {
            $jobLog{'exitStatus'} = $1;
          # Job used resources
          } elsif ($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=(\d+:\d+:\d+)/) {
            $jobLog{'cput'} = $1;
            $jobLog{'mem'} = $2;
            $jobLog{'vmem'} = $3;
            $jobLog{'walltime'} = $4;
            # Compute duration in seconds from walltime hours, minutes, seconds
            $jobLog{'duration'} = timeToSeconds($4);
          # Job end date
          } elsif ($jobLine =~ /^End PBS Epilogue (.*) (\d+)$/) {
            $jobLog{'endDate'} = $1;
            $jobLog{'endSecondsSinceEpoch'} = $2;
          }
        }
        close(CLUSTER_JOB_LOG_FILE);
      }
      push (@$rAoH_jobLogList, \%jobLog);
    }
  }
  close(JOB_LOG_LIST_FILE);
}

# Print out a log report in simple text format
sub getLogTextReport {

  my $jobLogListPath = shift;
  my @AoH_jobLogList;

  # Parse the job log files and get an array of hashes of those logs
  readJobLogListFile($jobLogListPath, \@AoH_jobLogList);

  my $logTextReport = "";

  $logTextReport .= "# Number of jobs: " . ($#AoH_jobLogList + 1) . "\n#\n";

  # Retrieve first job start date, last job end date, shortest job, longest job
  my $firstStartSecondsSinceEpoch;
  my $lastEndSecondsSinceEpoch;
  my $shortestJob;
  my $longestJob;

  for my $jobLog (@AoH_jobLogList) {
    if (exists $jobLog->{'startSecondsSinceEpoch'} and (not defined $firstStartSecondsSinceEpoch or $firstStartSecondsSinceEpoch > $jobLog->{'startSecondsSinceEpoch'})) {
      $firstStartSecondsSinceEpoch = $jobLog->{'startSecondsSinceEpoch'};
    }
    if (exists $jobLog->{'endSecondsSinceEpoch'} and (not defined $lastEndSecondsSinceEpoch or $lastEndSecondsSinceEpoch < $jobLog->{'endSecondsSinceEpoch'})) {
      $lastEndSecondsSinceEpoch = $jobLog->{'endSecondsSinceEpoch'};
    }
    if (exists $jobLog->{'duration'}) {
      if (not defined $shortestJob or $shortestJob->{'duration'} > $jobLog->{'duration'}) {
        $shortestJob = $jobLog;
      }
      if (not defined $longestJob or $longestJob->{'duration'} < $jobLog->{'duration'}) {
        $longestJob = $jobLog;
      }
    }
  }

  # Print out execution time
  my $executionTime = (defined $firstStartSecondsSinceEpoch and defined $lastEndSecondsSinceEpoch) ? formatDuration($lastEndSecondsSinceEpoch - $firstStartSecondsSinceEpoch) : "N/A";
  my $startDate = defined $firstStartSecondsSinceEpoch ? strftime('%FT%T', localtime($firstStartSecondsSinceEpoch)) : "N/A";
  my $endDate = defined $lastEndSecondsSinceEpoch ? strftime('%FT%T', localtime($lastEndSecondsSinceEpoch)) : "N/A";
  $logTextReport .= "# Execution time: $startDate - $endDate ($executionTime)\n#\n";

  # Print out shortest and longest jobs
  $logTextReport .= "# Shortest job: " . (defined $shortestJob ? $shortestJob->{'jobName'} . " (" . formatDuration($shortestJob->{'duration'}) . ")" : "N/A") . "\n";
  $logTextReport .= "# Longest job: " . (defined $longestJob ? $longestJob->{'jobName'} . " (" . formatDuration($longestJob->{'duration'}) . ")" : "N/A") . "\n#\n";

  $logTextReport .= join("\t", (
    "#JOB_ID",
    "JOB_FULL_ID",
    "JOB_NAME",
    "JOB_DEPENDENCIES",
    "JOB_EXIT_CODE",
    "CMD_EXIT_CODE",
    "REAL_TIME",
    "START_DATE",
    "END_DATE",
    "CPU_TIME",
    "CPU_REAL_TIME_RATIO",
    "PHYSICAL_MEM",
    "VIRTUAL_MEM",
    "EXTRA_VIRTUAL_MEM_PCT",
    "PATH"
  )) . "\n";

  for my $jobLog (@AoH_jobLogList) {
    $logTextReport .= join("\t", (
      exists $jobLog->{'jobId'} ? $jobLog->{'jobId'} : "N/A",
      exists $jobLog->{'jobFullId'} ? $jobLog->{'jobFullId'} : "N/A",
      exists $jobLog->{'jobName'} ? $jobLog->{'jobName'} : "N/A",
      exists $jobLog->{'jobDependencies'} ? $jobLog->{'jobDependencies'} : "N/A",
      exists $jobLog->{'exitStatus'} ? $jobLog->{'exitStatus'} : "N/A",
      exists $jobLog->{'MUGQICexitStatus'} ? $jobLog->{'MUGQICexitStatus'} : "N/A",
      exists $jobLog->{'walltime'} ? $jobLog->{'walltime'} . " (" . formatDuration($jobLog->{'duration'}) . ")" : "N/A",
      exists $jobLog->{'startSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'startSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'endSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'endSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'cput'} ? $jobLog->{'cput'} . " (" . formatDuration(timeToSeconds($jobLog->{'cput'})) . ")" : "N/A",
      (exists $jobLog->{'walltime'} and exists $jobLog->{'cput'} and timeToSeconds($jobLog->{'walltime'}) != 0) ? sprintf("%.2f", timeToSeconds($jobLog->{'cput'}) / timeToSeconds($jobLog->{'walltime'})) : "N/A",
      exists $jobLog->{'mem'} ? sprintf("%.1f", kiBToGiB($jobLog->{'mem'})) . " GiB" : "N/A",
      exists $jobLog->{'vmem'} ? sprintf("%.1f", kiBToGiB($jobLog->{'vmem'})) . " GiB" : "N/A",
      (exists $jobLog->{'vmem'} and exists $jobLog->{'mem'} and kiBToGiB($jobLog->{'mem'}) != 0) ? sprintf("%.1f", (kiBToGiB($jobLog->{'vmem'}) / kiBToGiB($jobLog->{'mem'}) - 1) * 100) . " %" : "N/A",
      exists $jobLog->{'path'} ? $jobLog->{'path'} : "N/A"
    )) . "\n";
  }

  return $logTextReport;
}

# Return seconds from time given in hh:mm:ss
sub timeToSeconds {
  my $time =  shift;

  if ($time =~ /^(\d+):(\d\d):(\d\d)/) {
    return $1 * 60 ** 2 + $2 * 60 + $3;
  } else {
    return "N/A";
  }
}

# Return duration given in seconds into human readable format
sub formatDuration {
  my $seconds =  shift;

  # Less than 1 minute
  if ($seconds < 60) {
    return $seconds . " s";
  }
  # Less than 1 hour
  elsif ($seconds < (60 * 60)) {
    return int($seconds / 60 ) . " min " . formatDuration($seconds % 60);
  }
  # Less than 1 day
  elsif ($seconds < (60 * 60 * 24)) {
    return int($seconds / (60 * 60)) . " h " . formatDuration($seconds % (60 * 60));
  }
  # 1 day or more
  else {
    return int($seconds / (60 * 60 * 24)) . " d " . formatDuration($seconds % (60 * 60 * 24));
  }
}

sub kiBToGiB {
  my $size = shift;

  if ($size =~ /^(\d+)kb/) {
    return $1 / (1024 ** 2);
  } else {
    return "N/A";
  }
}
        
1;
