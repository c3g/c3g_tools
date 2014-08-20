#!/usr/bin/env perl

use strict;
use POSIX qw( strftime );
use Net::SMTP;
#use Time::localtime;
#use File::stat;
use Getopt::Long;
use Time::Piece;

my $version = "0.1";

&main();

sub main {
  my @dirs;
  my @to;
  my $verbose;
  my $fromEmail;
  my $server;
  my $smtpHost;
  my $result = GetOptions ( "from=s" => \$fromEmail,
                            "to=s" => \@to,
                            "dir=s" => \@dirs,
                            "server=s" => \$server,
                            "smtpHost=s" => \$smtpHost,
                            "verbose" => \$verbose);

  if(!defined($server) || length($server) == 0) {
    die ("You need to set server");
  }

  if(!defined($smtpHost) || length($smtpHost) == 0) {
    $smtpHost = 'mailserver.genome.mcgill.ca';
  }

  my %dirInfo;
  my %user2Dirs;
  my %userTotal;
  for my $rootDir (@dirs) {
    if(! -d $rootDir) {
      warn "'".$rootDir."' is not a directory. Ignored\n";
    }

    opendir(ROOT_DIR, $rootDir) or die "Couldn't open directory ".$rootDir."\n";
    my @dirsToTest =  grep { /^[^\.]/ && -d "$rootDir/$_" } readdir(ROOT_DIR);
    closedir(ROOT_DIR);

    for my $relDirToTest (@dirsToTest) {
      my $dirToTest = $rootDir.'/'.$relDirToTest;
      my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($dirToTest);
      my $uname = getpwuid($uid);
      my $du = `du -bs $dirToTest`;
      my @duInfos = split('\t', $du);
      my $timestamp       = localtime($mtime)->dmy;

      if(!defined($uname) || length($uname) == 0) {
        $uname == 'UNKNOWN';
      }

      $dirInfo{$dirToTest} = {'dir' => $dirToTest, 'uname' => $uname, 'du' => $duInfos[0], 'mtime' => $timestamp};
      if(!defined($user2Dirs{$uname})) {
        $user2Dirs{$uname} = [];
      }
      push(@{$user2Dirs{$uname}}, $dirInfo{$dirToTest});
      $userTotal{$uname} += $dirInfo{$dirToTest}->{'du'};
    }
  }

  my @userPerSize = sort { $userTotal{$b} <=> $userTotal{$a} } keys %userTotal;
  my @dirPerSize = sort { $dirInfo{$b}->{'du'} <=> $dirInfo{$a}->{'du'} } keys %dirInfo;

  my $output;
  $output .= 'Users:'."\n";
  $output .= '======'."\n";
  for my $user (@userPerSize) {
    $output .= sprintf("%-10s%-10s\n", $user, nice_size($userTotal{$user}));
  }
  $output .= "\n";
  $output .= 'Dir Per Users:'."\n";
  $output .= '=============='."\n";
  for my $user (@userPerSize) {
    $output .= $user."\n";
    $output .= '========='."\n";

    my @sortedUserDir = sort { $b->{'du'} <=> $a->{'du'} } @{$user2Dirs{$user}};
    for my $rH_userDir (@sortedUserDir) {
      $output .= sprintf("%-10s%-11s%-100s\n", nice_size($rH_userDir->{'du'}), $rH_userDir->{'mtime'}, $rH_userDir->{'dir'});
    }
    $output .= "\n";
  }
  $output .= "\n";
 
  $output .= "\n";
  $output .= 'Directories:'."\n";
  $output .= '============'."\n";
  for my $dir (@dirPerSize) {
    $output .= sprintf("%-10s%-10s%-11s%-100s\n", nice_size($dirInfo{$dir}->{'du'}), $dirInfo{$dir}->{'uname'}, $dirInfo{$dir}->{'mtime'}, $dir);
  }
  $output .= "\n";

  print $output;
  sendEmail($fromEmail, \@to, $server, $smtpHost, $output);
}

sub nice_size {
  my $size = shift;
  my $i = 0;

  while ($size > 1000) {
    $size = $size / 1000;
    $i++;
  }

  my @sizes=('B','KB','MB','GB','TB','PB');
  return sprintf("%.1f$sizes[$i]", $size);
}

sub sendEmail {
  my $fromEmail = shift;
  my $rA_to = shift;
  my $server = shift;
  my $smtpHost = shift;
  my $message = shift;

  my $smtp = Net::SMTP->new($smtpHost);
  if(!defined($smtp) || !($smtp)) {
    print STDERR "SMTP ERROR: Unable to open smtp session.\n";
    return 1;
  }

  if (! ($smtp->mail($fromEmail) ) ) {
    print STDERR "SMTP ERROR: Cannot set from: ".$fromEmail.".\n";
    return 2;
  }

  if (! ($smtp->to(@$rA_to) ) ) {
    print STDERR "SMTP ERROR: Cannot set to: ".join(',', @$rA_to).".\n";
    return 3;
  }

  my $tos = join(',',@$rA_to);

	my $subject = 'Subject: MUGQIC '.$server.' Space Usage Breakdown'."\n";

  $smtp->data();
  $smtp->datasend("To: ".$tos."\n");
  $smtp->datasend($subject);
  $smtp->datasend("\n");
  $smtp->datasend($message);
  $smtp->dataend();

  $smtp->quit;

  return 0;
}

