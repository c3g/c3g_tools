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
use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
getLogReport.pl

USAGE:
getLogReport.pl [options] <joblist file>

OPTIONS:
--memtime : Flag if memtime output if present in your jobs output files.
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Joel Fillon - joel.fillon@mcgill.ca

ENDHERE

# Check and assign the parameters
my($memtime, $help);
GetOptions( 
  "memtime" => \$memtime,
  "help"    => \$help
);
if ($help) { print $usage; exit; }

#my $jobLogList = unshift @ARGV;
my $jobLogList = $ARGV[0];

print getLogTextReport($jobLogList);

# SUB
#-----------------------

# Return an array of hashes of log values for each job's log file created by the cluster
sub readJobLogListFile {
  my $jobLogListPath = shift;
  # Return array of hash of log values
  my $rAoH_jobLogList = shift;
  my $rAoHoH_jobLogMemtime = shift;

  # Read the global list of log files
  open(JOB_LOG_LIST_FILE, $jobLogListPath) or die "Cannot open $jobLogListPath\n";
  while(my $line = <JOB_LOG_LIST_FILE>) {
    chomp($line);

    # Retrieve each job's log file path
    my ($jobId, $jobName, $jobDependencies, $clusterJobLogPath) = split(/\t/, $line);

    if (defined $jobName and defined $clusterJobLogPath) {
      my %jobLog;
      my %memtime;

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
		my $i = 0;
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
          # Username
          } elsif ($jobLine =~ /^Username:\s+(\S+)/) {
            $jobLog{'username'} = $1;
          # Group
          } elsif ($jobLine =~ /^Group:\s+(\S+)/) {
            $jobLog{'group'} = $1;
          # Session
          } elsif ($jobLine =~ /^Session:\s+(\S+)/) {
            $jobLog{'session'} = $1;
          # Limits
          } elsif ($jobLine =~ /^Limits:\s+(\S+)/) {
            $jobLog{'limits'} = $1;
          # Job used resources
          } elsif ($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=(\d+:\d+:\d+)/) {
            $jobLog{'cput'} = $1;
            $jobLog{'mem'} = $2;
            $jobLog{'vmem'} = $3;
            $jobLog{'walltime'} = $4;
            # Compute duration in seconds from walltime hours, minutes, seconds
            $jobLog{'duration'} = timeToSeconds($4);
          # Queue
          } elsif ($jobLine =~ /^Queue:\s+(\S+)/) {
            $jobLog{'queue'} = $1;
          # Account
          } elsif ($jobLine =~ /^Account:\s+(\S+)/) {
            $jobLog{'account'} = $1;
          # Job exit status (should be the same as MUGQIC exit status unless MUGQIC exit status is skipped)
          } elsif ($jobLine =~ /^Exit_status:\s+(\d+)/) {
            $jobLog{'exitStatus'} = $1;
          # Nodes
          } elsif ($jobLine =~ /^Nodes:\s+(\S+)/) {
            $jobLog{'nodes'} = $1;
          # Job end date
          } elsif ($jobLine =~ /^End PBS Epilogue (.*) (\d+)$/) {
            $jobLog{'endDate'} = $1;
            $jobLog{'endSecondsSinceEpoch'} = $2;
          } 
          
          # Memtime if present in options (-m).
          if($memtime){
			if ($jobLine =~ m/(\d+\.\d+) user, (\d+\.\d+) system, (\d+\.\d+) elapsed -- Max VSize = (\d+)KB, Max RSS = (\d+)KB/) {
              $memtime{$i}{'memTime_user'} = $1;
           	  $memtime{$i}{'memTime_system'} = $2;
           	  $memtime{$i}{'memTime_elapsed'} = $3;
           	  $memtime{$i}{'memTime_MaxVSize'} = $4;
           	  $memtime{$i}{'memTime_MaxRSS'} = $5;
              #print STDERR "[DEBUG]\t\$i:".$i."\n";
              $i++;
              #print STDERR "[DEBUG]\t1:\t".$1."\n"."2:\t".$2."\n"."3:\t".$3."\n"."3:\t".$3."\n"."3:\t".$3."\n"."4:\t".$4."\n"."5:\t".$5."\n";
            }
          }
        }
        close(CLUSTER_JOB_LOG_FILE);
      }
      push (@$rAoH_jobLogList, \%jobLog);
      push (@$rAoHoH_jobLogMemtime, \%memtime);
    }
  }
  close(JOB_LOG_LIST_FILE);
}

# Print out a log report in simple text format
sub getLogTextReport {

  my $jobLogListPath = shift;
  my @AoH_jobLogList;
  my @AoHoH_jobLogMemtime;

  # Parse the job log files and get an array of hashes of those logs
  readJobLogListFile($jobLogListPath, \@AoH_jobLogList, \@AoHoH_jobLogMemtime);

  my $logTextReport = "";

  $logTextReport .= "# Number of jobs: " . ($#AoH_jobLogList + 1) . "\n#\n";

  # Retrieve first job start date, last job end date, shortest/longest jobs, lowest/highest memory jobs
  my $firstStartSecondsSinceEpoch;
  my $lastEndSecondsSinceEpoch;
  my $shortestJob;
  my $longestJob;
  my $lowestMemoryJob;
  my $highestMemoryJob;

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
    if (exists $jobLog->{'mem'}) {
      if (not defined $lowestMemoryJob or kiBToNum($lowestMemoryJob->{'mem'}) > kiBToNum($jobLog->{'mem'})) {
        $lowestMemoryJob = $jobLog;
      }
      if (not defined $highestMemoryJob or kiBToNum($highestMemoryJob->{'mem'}) < kiBToNum($jobLog->{'mem'})) {
        $highestMemoryJob = $jobLog;
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
  $logTextReport .= "# Longest job: " . (defined $longestJob ? $longestJob->{'jobName'} . " (" . formatDuration($longestJob->{'duration'}) . ")" : "N/A") . "\n";
  $logTextReport .= "#\n";

  # Print out lowest and highest memory jobs
  $logTextReport .= "# Lowest memory job: " . (defined $lowestMemoryJob ? $lowestMemoryJob->{'jobName'} . " (" . sprintf("%.2f", kiBToGiB($lowestMemoryJob->{'mem'})) . " GiB" . ")" : "N/A") . "\n";
  $logTextReport .= "# Highest memory job: " . (defined $highestMemoryJob ? $highestMemoryJob->{'jobName'} . " (" . sprintf("%.2f", kiBToGiB($highestMemoryJob->{'mem'})) . " GiB" . ")" : "N/A") . "\n";
  $logTextReport .= "#\n";

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
    "LIMITS",
    "QUEUE",
    "USERNAME",
    "GROUP",
    "SESSION",
    "ACCOUNT",
    "NODES",
    "PATH"
  ));

  if($memtime){
    $logTextReport .= join("\t", (
      "MEMTIME_USER (hrs)",
      "MEMTIME_SYSTEM (hrs)",
      "MEMTIME_ELAPSED (Gb)",
      "MEMTIME_MAXVSIZE (Gb)",
      "MEMTIME_MAXRSS (Gb)"
    )) . "\n";
  }
  
  for my $jobLog (@AoH_jobLogList) {
    my $memtimeLog = shift(@AoHoH_jobLogMemtime) if($memtime);

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
      exists $jobLog->{'mem'} ? sprintf("%.2f", kiBToGiB($jobLog->{'mem'})) . " GiB" : "N/A",
      exists $jobLog->{'vmem'} ? sprintf("%.2f", kiBToGiB($jobLog->{'vmem'})) . " GiB" : "N/A",
      (exists $jobLog->{'vmem'} and exists $jobLog->{'mem'} and kiBToGiB($jobLog->{'mem'}) != 0) ? sprintf("%.1f", (kiBToGiB($jobLog->{'vmem'}) / kiBToGiB($jobLog->{'mem'}) - 1) * 100) . " %" : "N/A",
      exists $jobLog->{'limits'} ? $jobLog->{'limits'} : "N/A",
      exists $jobLog->{'queue'} ? $jobLog->{'queue'} : "N/A",
      exists $jobLog->{'username'} ? $jobLog->{'username'} : "N/A",
      exists $jobLog->{'group'} ? $jobLog->{'group'} : "N/A",
      exists $jobLog->{'session'} ? $jobLog->{'session'} : "N/A",
      exists $jobLog->{'account'} ? $jobLog->{'account'} : "N/A",
      exists $jobLog->{'nodes'} ? $jobLog->{'nodes'} : "N/A",
      exists $jobLog->{'path'} ? $jobLog->{'path'} : "N/A"
    ));

	if($memtime){
		$logTextReport .= "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n";
	}else{
		$logTextReport .= "\n";
	}
  
    # memtime log on a separate line.
    if($memtime){
      foreach my $key (sort{$a <=> $b} keys %$memtimeLog) {
		#print STDERR "[DEBUG]\t".$key."\n";
		#print STDERR $memtimeLog->{$key}{'memTime_user'}."\n";
	
        $logTextReport .= join("\t", (
          "N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A",
          exists $memtimeLog->{$key}{'memTime_user'} ?  secondsToHours($memtimeLog->{$key}{'memTime_user'}) : "N/A",
          exists $memtimeLog->{$key}{'memTime_system'} ? secondsToHours($memtimeLog->{$key}{'memTime_system'}) : "N/A",
          exists $memtimeLog->{$key}{'memTime_elapsed'} ? secondsToHours($memtimeLog->{$key}{'memTime_elapsed'}) : "N/A",
          exists $memtimeLog->{$key}{'memTime_MaxVSize'} ? sprintf("%.6f", kiBToGiB($memtimeLog->{$key}{'memTime_MaxVSize'} )) : "N/A",
          exists $memtimeLog->{$key}{'memTime_MaxRSS'} ? sprintf("%.6f", kiBToGiB($memtimeLog->{$key}{'memTime_MaxRSS'} )) : "N/A"
        )) . "\n";
      }
    }
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

# Return hours from time given seconds in 00.00 (say 10 sec and 10/100 secs...)
sub secondsToHours {
  my $time =  shift;

  if ($time =~ /^(\d+)\.(\d+)/) {
    return sprintf("%.2f", ($1 / 60 /  60));
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

sub kiBToNum {
  my $size = shift;

  if ($size =~ /^(\d+)kb/i) {
    return $1;
  } else {
    return "N/A";
  }
}

sub kiBToGiB {
  my $size = shift;

  if ($size =~ /^(\d+)kb/) {
    return $1 / (1024 ** 2);
  } elsif ($size =~ /^(\d+)/){
    return $1 / (1024 ** 2);
  }else {
    return "N/A";
  }
}
        
1;
