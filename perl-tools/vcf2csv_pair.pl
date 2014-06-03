#!/usr/bin/perl

use Time::localtime;
use List::Util qw(sum);
use strict;

my $EXPERIMENT_TYPE = $ENV{EXPERIMENT_TYPE};

&main();

my @chromosomes;
my $currentChr;
sub parseChr {
  my $chrFile = shift;
  open(CHR, $chrFile) or die "Cannot open $chrFile\n";
  while(<CHR>) {
    chomp($_);
    #my @chr = split(/\t/, $_);
    #my $chrom = $chr[1];
    my $chrom = $_;
    $chrom =~ s/chr//g;
    push(@chromosomes, $chrom);
  }
  close(CHR);
  $currentChr=shift(@chromosomes);
  print STDERR "Step 2: Chromosome file loaded\n";
}

sub initVcfFiles {
  my $rA_vcfFiles = shift;  
  my $rH_normals = shift;
  my $rH_tumors = shift;

  my %fileDetails;
  my %samples;
  for my $vcfFile (@{$rA_vcfFiles}) {
    #my ($sampleName) = $vcfFile =~ /^([^.\/]+)/;
    my ($sampleName) = $vcfFile =~ /^.*\/(.*)\.merged.*/;
    my $fh;
    open($fh, "<", $vcfFile) or die "Can't open $vcfFile\n";
    my $line;
    while($line = <$fh>) {
      if($line =~ /^#CHR/) {
        last;
      }
    }
    chomp($line);
    my @values = split(/\t/,$line);
    if(@values < 11) {
      warn "File $vcfFile doesn't have enough columns\n";
      close($fh);
      next;
    }
    my $normalSample = $values[9];
    my $tumorSample  = $values[10];
    if( !defined( $rH_normals->{$normalSample} ) ) {
      warn "Normal sample names don't match: $sampleName\n";
#      next;
    }
    ;
    if(!defined( $rH_tumors->{$tumorSample} ) ) {
      warn "Tumor sample names don't match: $sampleName\n";
#      next;
    }

    $fileDetails{$sampleName} = {};
    $fileDetails{$sampleName}->{'fh'} = $fh;
    $line = <$fh>;
    chomp($line);
    $fileDetails{$sampleName}->{'currentLine'} = $line;
  }
  print STDERR "STEP 4: Check vcf files\n";
  return \%fileDetails;
}

sub nextSamples {
  my $rHoH_samplesFiles = shift;

  my $start=9999999999;

  my $chromosomeMatch=0;
  my %samplesToPrint;
  for my $sample (keys(%{$rHoH_samplesFiles})) {
    my $line = $rHoH_samplesFiles->{$sample}->{'currentLine'};
    my @lineValues = split(/\t/, $line);
    my $chromosome = $lineValues[0];
    my $position = $lineValues[1];

    if($chromosome eq $currentChr && $position <= $start){
      $chromosomeMatch=1;
      if($position < $start) {
        %samplesToPrint = ();
      }
      $samplesToPrint{$sample} = \@lineValues;
      $start=$position;
    }
  }

  if($chromosomeMatch == 0) {
    $currentChr=shift(@chromosomes);
    if(!defined($currentChr)) {
      return undef;
    }
    return nextSamples($rHoH_samplesFiles);
  }
  else {
    for my $sample (keys(%samplesToPrint)) {
      my $fh = $rHoH_samplesFiles->{$sample}->{'fh'};
      my $line = <$fh>;
      if(!defined($line)) {
        close($fh);
        delete $rHoH_samplesFiles->{$sample};
      }
      else {
        chomp($line);
        $rHoH_samplesFiles->{$sample}->{'currentLine'} = $line;
      }
    }
    return \%samplesToPrint;
  }
}

sub main {
  my ($chromosomeFile, $normalTumorFile, @vcfFiles) = @ARGV;
  #my ($chromosomeFile, $kgFile, $normalTumorFile, @vcfFiles) = @ARGV;

  #my $kgFh;
  #open($kgFh, $kgFile) or die "Can't open $kgFile\n";
  #init1kgMAF($kgFh);

  parseChr($chromosomeFile);

  my ($rH_normals, $rH_tumors, $rHoA_alias) = parseSampleNormalTumor($normalTumorFile);

  my $rHoH_samplesFiles = initVcfFiles(\@vcfFiles, $rH_normals, $rH_tumors);
  my @samplesGiven = keys(%{$rHoH_samplesFiles});

  print "chromosome\tposition\tref_allele\talt_alleles\tAvg. Variant Quality\tmappability_flag\tgene_name"; #gene_length(aa)";
  print "\ttumor:normal_ratio";

  #print"gene_description\tgo_ids\tgo_terms\t";
  print "\t1000_genome_global_maf\timpact\teffect\teffect_type\tcodon_change\tamino_acid_change\tgene_type\ttranscript\tdb_snp_id\tCosmicID\t";
  print "InterPro Domain\tUniprot\tGERP++ neutral rate\tGERP++ RS score\t29way SiPhy score\tPolyphen2 prediction based on HumVar, 'D' ('porobably damaging') 'P' ('possibly damaging') and 'B' ('benign')\tSIFT score\t";
  print "Normal HomRef\tNormal HomAlt\tNormal Het\tTumor HomRef\tTumor HomAlt\tTumor Het\t";
  for my $sampleName (sort @samplesGiven) {
    print $sampleName." Variant quality\t";
    print $sampleName." CLR (phred pair GT likelyhood)\t";
    print $sampleName." Depth\t";
    print $sampleName."\t";
  }
  print "ucsc_url\tfull_effect_description\n"; 

  while(1) {
    my $rHoA_samplesToPrint = nextSamples($rHoH_samplesFiles);
    if(!defined($rHoA_samplesToPrint)) {
      last;
    }

    my %snps;
    my %sample_count = ();
    my %samples = ();
    
    for my $sample (keys(%{$rHoA_samplesToPrint})) {
      getSnpInformation($sample, \%snps, $rHoA_samplesToPrint, $rH_normals, $rH_tumors, $rHoA_alias, \%sample_count, \%samples);
      #getSnpInformation($sample, \%snps, $kgFh, $rHoA_samplesToPrint, $rH_normals, $rH_tumors, $rHoA_alias, \%sample_count, \%samples);
    }
    
    printSnps(\@samplesGiven, \%snps, $rH_normals, $rH_tumors, $rHoA_alias, \%sample_count, \%samples);
  }
}

sub getSnpInformation {
  my $sample = shift;
  my $rHoH_snps = shift;
  #my $kgFh = shift;
  my $rHoA_samplesToPrint = shift;
  my $rH_normals = shift;
  my $rH_tumors = shift;
  my $rHoA_alias = shift;
  my $rH_sampleCount = shift;
  my $rHoA_samples = shift;

  my $chromosome = $rHoA_samplesToPrint->{$sample}->[0];
  my $position = $rHoA_samplesToPrint->{$sample}->[1];
  my $knownSnpID = $rHoA_samplesToPrint->{$sample}->[2];
  my $refAllele = $rHoA_samplesToPrint->{$sample}->[3];
  my $altAlleles = $rHoA_samplesToPrint->{$sample}->[4];
  my $varQual = $rHoA_samplesToPrint->{$sample}->[5];
  my $genotypeQuality;
  my $mappabilityFlag;
  my $coverage;
  my $nmbRef;
  my $nmbAlt;

  my %snp;
  $snp{"chr"} = $chromosome;
  $snp{"position"} = $position;

  my @readCountsForSnp;
  my $fullEffect;
  my %effects;
  $effects{"HIGH"} = [];
  $effects{"MODERATE"} = [];
  $effects{"LOW"} = [];
  $effects{"OTHER"} = [];

  my @infos = split(/;/, $rHoA_samplesToPrint->{$sample}->[7]);
  my $clr = undef;
  for my $info (@infos) {
    if($info =~ /^EFF=(.*)/) {
      $fullEffect = $1;
      my @transcriptsEffect = split(/,/, $fullEffect);
      for my $transcriptEffect (@transcriptsEffect) {
        my %effectDetails;
        my ($effect, $details) = $transcriptEffect =~ /([^\(]+)\((.+)\)/;

        my @effInfos = split(/\|/, $details);
        $effectDetails{"effect"} = $effect;
        $effectDetails{"effectImpact"} = $effInfos[0];
        $effectDetails{"effectType"} = $effInfos[1];
        $effectDetails{"codonChange"} = $effInfos[2];
        $effectDetails{"aaChange"} = $effInfos[3];
	$effectDetails{"aaLength"} = $effInfos[4];   
        $effectDetails{"gene"} = $effInfos[5];
        $effectDetails{"bioType"} = $effInfos[6];
        $effectDetails{"coding"} = $effInfos[7];
        $effectDetails{"transcript"} = $effInfos[8];
        $effectDetails{"numExon"} = $effInfos[9];

        if($effectDetails{"effectImpact"} eq "HIGH") {
          push(@{$effects{"HIGH"}}, \%effectDetails);
        }
        elsif($effectDetails{"effectImpact"} eq "MODERATE") {
          push(@{$effects{"MODERATE"}}, \%effectDetails);
        }
        elsif($effectDetails{"effectImpact"} eq "LOW") {
          # not needed in snpEffect > 2.0.4_rc1
          if( $effectDetails{"effect"} eq "INTRAGENIC" ||
              $effectDetails{"effect"} eq "INTERGENIC" || 
              $effectDetails{"effect"} eq "INTERGENIC_CONSERVED" || 
              $effectDetails{"effect"} eq "UPSTREAM" || 
              $effectDetails{"effect"} eq "INTRON") {
            $effectDetails{"effectImpact"} = 'MODIFIER';
            push(@{$effects{"OTHER"}}, \%effectDetails);
          }
          else {
            push(@{$effects{"LOW"}}, \%effectDetails);
          }
        }
        else {
          push(@{$effects{"OTHER"}}, \%effectDetails);
        }
      }
    }
    elsif($info =~ /^DP4=(.*)/) {
      @readCountsForSnp = split(/,/,$1);
    }
    elsif($info =~ /^MIL=(.*)/) {
      $snp{"mappabilityFlag"} = $1;
    }
    #elsif($info =~ /^GENE_ANNOT=(.*)/) {
    #  my $geneDescriptions = $1;
      
    #  my @geneInfos = $geneDescriptions =~ m/(\([^|]+\|[^|]+\|[^|]+\|(?:\([^\)]+\),?)*\),?)/g;
    #  for my $geneInfo (@geneInfos) {
    #    my ($transcript, $geneName, $geneDesc, $goInfos) = $geneInfo =~ m/\(([^|]+)\|([^|]+)\|([^|]+)\|(.*)/g;
    #    if(defined($snp{"geneDescription"})) {
    #      $snp{"geneDescription"} = $snp{"geneDescription"} . ',';
    #    }
    #    $snp{"geneDescription"} .= $geneDesc;

    #    if(defined($snp{"goTerms"})) {
    #      $snp{"goTerms"} = $snp{"goTerms"} . ',';
    #      $snp{"goIds"} = $snp{"goIds"} . ',';
    #    }

    #    if(defined($goInfos) && length($goInfos) > 0) {
    #      my @goDetails = $goInfos =~ m/\(([^|]+)\|([^|]+)\|([^|]+)\),?/g;
          
    #      for(my $idx=0;$idx < @goDetails; $idx +=3) {
    #        $snp{"goIds"} .= '|'.$goDetails[$idx+0];
    #        $snp{"goTerms"} .= '|'.$goDetails[$idx+1];
    #      }
    #    }
    #  }
    #}
    elsif($info =~ /^CAF=\[[0-9.]+,([0-9.]+)/) {
      $snp{"GMAF"} = $1;
    }
    elsif($info =~ /^GMAF=([0-9.]+)/) {
      $snp{"GMAF"} = $1;
    }
    elsif($info =~ /^dbNSFP_GERP\+\+_NR=(.*)/) {
      $snp{"gerpNR"} = $1;
    }
    elsif($info =~ /^dbNSFP_GERP\+\+_RS=(.*)/) {
      $snp{"gerpRS"} = $1;
    }
    elsif($info =~ /^dbNSFP_29way_logOdds=(.*)/) {
      $snp{"29way"} = $1;
    }
    elsif($info =~ /^dbNSFP_Polyphen2_HVAR_pred=(.*)/) {
      $snp{"polyphenHvar"} = $1;
    }
    elsif($info =~ /^dbNSFP_Uniprot_acc=(.*)/) {
      $snp{"uniprot"} = $1;
    }
    elsif($info =~ /^dbNSFP_Interpro_domain=(.*)/) {
      $snp{"interPro"} = $1;
    }
    elsif($info =~ /^dbNSFP_SIFT_score=(.*)/) {
      $snp{"sift_score"} = $1;
    }    
    elsif($info =~ /^COSMIC_ID=(.*)/) {
      $snp{"cosmic"} = $1;
    }
    elsif($info =~ /^CLR=(.*)/) {
      $clr = $1;
    }    
  } # Variant INFO column

#      my @sampleFormat = split(/:/, $rHoA_samplesToPrint->{$sample}->[8]);
#      my $genotypeQualityIdx=-1;
#      for(my $idx=0 ; $idx < @sampleFormat; $idx++) {
#        if($sampleFormat[$idx] eq "GQ") {
#          $genotypeQualityIdx=$idx;
#        }
#      }
#      my @sampleInfos = split(/:/, $rHoA_samplesToPrint->{$sample}->[9]);
#      my $genotypeQuality = $sampleInfos[$genotypeQualityIdx];

  my $rH_defaultEffect;
  if(@{$effects{"HIGH"}} > 0) {
    $rH_defaultEffect = $effects{"HIGH"}->[0];
  }
  elsif(@{$effects{"MODERATE"}} > 0) {
    $rH_defaultEffect = $effects{"MODERATE"}->[0];
  }
  elsif(@{$effects{"LOW"}} > 0) {
    $rH_defaultEffect = $effects{"LOW"}->[0];
  }
  else {
    $rH_defaultEffect = $effects{"OTHER"}->[0];
  }

  $snp{"gene"} = $rH_defaultEffect->{"gene"};
  $snp{"impact"} = $rH_defaultEffect->{"effectImpact"};
  $snp{"effect"} = $rH_defaultEffect->{"effect"};
  $snp{"effectType"} = $rH_defaultEffect->{"effectType"};
  $snp{"codonChange"} = $rH_defaultEffect->{"codonChange"};
  $snp{"aaChange"} = $rH_defaultEffect->{"aaChange"};
  $snp{"aaLength"} = $rH_defaultEffect->{"aaLength"};
  $snp{"bioType"} = $rH_defaultEffect->{"bioType"};
  $snp{"transcript"} = $rH_defaultEffect->{"transcript"};
  $snp{"refAllele"} = $refAllele;
  $snp{"altAllele"} = {};
  $snp{"altAllele"}->{$altAlleles} = 1;
  $snp{"dbSnpID"} = $knownSnpID;
  $snp{"fullEffect"} = $fullEffect;
  #my $rH_1kgMAFs = get1kgMAF($snp{"chr"},$snp{"position"});
  #    if(defined $rH_1kgMAFs) {
  #      $snp{"1kg_global_maf"} = $rH_1kgMAFs->{"1kg_global_maf"};
  #      $snp{"1kg_asian_maf"} = $rH_1kgMAFs->{"1kg_asian_maf"};
  #      $snp{"1kg_american_maf"} = $rH_1kgMAFs->{"1kg_american_maf"};
  #      $snp{"1kg_african_maf"} = $rH_1kgMAFs->{"1kg_african_maf"};
  #      $snp{"1kg_european_maf"} = $rH_1kgMAFs->{"1kg_european_maf"};
  #    }

  my $key = $chromosome.'-'.$position;
  
  $rH_sampleCount->{"$key"}++;
  push( @{ $rHoA_samples->{"$key"} }, $sample );

  if($infos[0] =~ /INDEL/) {
    $key .= "-INDEL";
  }
  if(!defined($rHoH_snps->{$key})) {
    $snp{"samples"} = {};
    $snp{"tumors"} = {};
    $snp{"normals"} = {};
    $snp{"avgVarQual"} = [$varQual];
    $rHoH_snps->{$key} = \%snp;
  }
  else {
    push(@{$snp{"avgVarQual"}}, $varQual);
    if(!($infos[0] =~ /INDEL/)) {
      if($rHoH_snps->{"$key"}->{"refAllele"} ne $snp{"refAllele"}) {
        warn "Ref alleles don't match: $key, $sample\n";
      }

      my @alleles = keys(%{ $snp{"altAllele"} });
      my $altAllele = $alleles[0];
      if(!defined( $rHoH_snps->{"$key"}->{"altAllele"}->{"$altAllele"}) ) {
	#warn "Info: ".$infos[0]."\n";
	#warn "Alt alleles don't match: $key, $sampleName\n";
        $rHoH_snps->{"$key"}->{"altAllele"}->{"$altAllele"} = 1;
      }
    }
  }

  my @normalSampleInfos = split(/:/, $rHoA_samplesToPrint->{$sample}->[9]);
  my @tumorSampleInfos = split(/:/, $rHoA_samplesToPrint->{$sample}->[10]);
  my $normSample = $rHoA_alias->{$sample}->[0];
  my $tumorSample = $rHoA_alias->{$sample}->[1];

  # test DP for tumor depth
#  if($tumorSampleInfos[2] > 9) {
    $rHoH_snps->{"$key"}->{"enoughDepth"} = 1;
#  }

  $rHoH_snps->{"$key"}->{"samples"}->{$normSample} = $rHoA_samplesToPrint->{$sample}->[9];
  $rHoH_snps->{"$key"}->{"samples"}->{$tumorSample} = $rHoA_samplesToPrint->{$sample}->[10];
  $rHoH_snps->{"$key"}->{"samples"}->{$normSample."DP"} = $normalSampleInfos[2];
  $rHoH_snps->{"$key"}->{"samples"}->{$tumorSample."DP"} = $tumorSampleInfos[2];
  $rHoH_snps->{"$key"}->{"samples"}->{$sample."CLR"} = $clr;
  $rHoH_snps->{"$key"}->{"samples"}->{$sample."QUAL"} = $varQual;

  if($normalSampleInfos[0] ne '0/0') {
    $rHoH_snps->{"$key"}->{"normals"}->{$normSample.'#'.$sample} = 1;
    if(substr($normalSampleInfos[0],0,1) ne substr($normalSampleInfos[0],2,1)) {
      $rHoH_snps->{"$key"}->{"normalsHet"} = $rHoH_snps->{"$key"}->{"normalsHet"} + 1;
    } 
    else {
      $rHoH_snps->{"$key"}->{"normalsHomAlt"} = $rHoH_snps->{"$key"}->{"normalsHomAlt"} + 1;
    }
  }
  else {
	$rHoH_snps->{"$key"}->{"normalsHomRef"} = $rHoH_snps->{"$key"}->{"normalsHomRef"} + 1;
  }

  if($tumorSampleInfos[0] ne '0/0') {
    $rHoH_snps->{"$key"}->{"tumors"}->{$tumorSample.'#'.$sample} = 1;
    if(substr($tumorSampleInfos[0],0,1) ne substr($tumorSampleInfos[0],2,1)) {
      $rHoH_snps->{"$key"}->{"tumorsHet"} = $rHoH_snps->{"$key"}->{"tumorsHet"} + 1;
    } 
    else {
      $rHoH_snps->{"$key"}->{"tumorsHomAlt"} = $rHoH_snps->{"$key"}->{"tumorsHomAlt"} + 1;
    }
  }
  else {
	$rHoH_snps->{"$key"}->{"tumorsHomRef"} = $rHoH_snps->{"$key"}->{"tumorsHomRef"} + 1;
  }

  my @alleles = keys(%{ $snp{"altAllele"} });
  my @trueAlleles = split(/,/, $alleles[0]);

  my $newGenotype = $tumorSampleInfos[0];
  $newGenotype =~ s/0/$snp{"refAllele"}/g;
  for(my $idx=0; $idx < @trueAlleles; $idx++) {
    my $off = $idx+1;
    $newGenotype =~ s/$off/$trueAlleles[$idx]/g;
  }
  $rHoA_samplesToPrint->{$sample}->[10] =~ s/$tumorSampleInfos[0]/$newGenotype/g;

  $newGenotype = $normalSampleInfos[0];
  $newGenotype =~ s/0/$snp{"refAllele"}/g;
  for(my $idx=0; $idx < @trueAlleles; $idx++) {
    my $off = $idx+1;
    $newGenotype =~ s/$off/$trueAlleles[$idx]/g;
  }
  $rHoA_samplesToPrint->{$sample}->[9] =~ s/$normalSampleInfos[0]/$newGenotype/g;

  #add1kgMAF($kgFh, $rHoH_snps, $key);
}

sub printSnps {
  my $rA_samplesGiven = shift;
  my $rHoH_snps = shift;
  my $rH_normals = shift;
  my $rH_tumors = shift;
  my $rHoA_alias = shift;
  my $rH_sampleCount = shift;
  my $rHoA_samples = shift;

  for my $rH_snp (values(%{$rHoH_snps})) {
    # skip low tumor depth variants
    if(!defined($rH_snp->{"enoughDepth"}) || $rH_snp->{"enoughDepth"} != 1) {
      next;
    }

    my $key = $rH_snp->{"chr"}."-".$rH_snp->{"position"};

    print $rH_snp->{"chr"};
    print "\t".$rH_snp->{"position"};
    print "\t".$rH_snp->{"refAllele"};
    print "\t".join('/', keys(%{$rH_snp->{"altAllele"}}));
    print "\t".(sum(@{$rH_snp->{"avgVarQual"}})/@{$rH_snp->{"avgVarQual"}});
    print "\t".$rH_snp->{"mappabilityFlag"};
    print "\t".$rH_snp->{"gene"};
    #print "\t".$rH_snp->{"aaLength"};
    print "\t\"".keys(%{$rH_snp->{"tumors"}}).":".keys(%{$rH_snp->{"normals"}})."\"";
    print "\t".$rH_snp->{"GMAF"};
    #print "\t".$rH_snp->{"geneDescription"};
    #print "\t".$rH_snp->{"goIds"};
    #print "\t".$rH_snp->{"goTerms"};
    print "\t".$rH_snp->{"impact"};
    print "\t".$rH_snp->{"effect"};
    print "\t".$rH_snp->{"effectType"};
    print "\t".$rH_snp->{"codonChange"};
    print "\t".$rH_snp->{"aaChange"};
    print "\t".$rH_snp->{"bioType"};
    print "\t".$rH_snp->{"transcript"};
    print "\t".$rH_snp->{"dbSnpID"};
    #print "\tASN:".$rH_snp->{"1kg_asian_maf"}.' AMR:'.$rH_snp->{"1kg_american_maf"}.' AFR:'.$rH_snp->{"1kg_african_maf"}.' EUR:'.$rH_snp->{"1kg_european_maf"};
    print "\t".$rH_snp->{"cosmic"};
    print "\t".$rH_snp->{"interPro"};
    print "\t".$rH_snp->{"uniprot"};
    print "\t".$rH_snp->{"gerpNR"};
    print "\t".$rH_snp->{"gerpRS"};
    print "\t".$rH_snp->{"29way"};
    print "\t".$rH_snp->{"polyphenHvar"};
    print "\t".$rH_snp->{"sift_score"};
    print "\t".(0+$rH_snp->{"normalsHomRef"});
    print "\t".(0+$rH_snp->{"normalsHomAlt"});
    print "\t".(0+$rH_snp->{"normalsHet"});
    print "\t".(0+$rH_snp->{"tumorsHomRef"});
    print "\t".(0+$rH_snp->{"tumorsHomAlt"});
    print "\t".(0+$rH_snp->{"tumorsHet"});
    print "\t";
   	for my $sampleName (sort @{$rA_samplesGiven}) {
   	  my $normSampleName = $rHoA_alias->{$sampleName}[0];
      my $tumorSampleName = $rHoA_alias->{$sampleName}[1];
   	  
   	  if(defined($rH_snp->{"samples"}->{$tumorSampleName})) {
   	    print $rH_snp->{"samples"}->{$sampleName."QUAL"}."\t";
   	    print $rH_snp->{"samples"}->{$sampleName."CLR"}."\t";
   	    print $rH_snp->{"samples"}->{$normSampleName."DP"}." , ".$rH_snp->{"samples"}->{$tumorSampleName."DP"}."\t";
   	    print $rH_snp->{"samples"}->{$normSampleName}." , ".$rH_snp->{"samples"}->{$tumorSampleName}."\t";
   	  }
      else {
        print "\t";
        print "\t";
        print "\t";
        print "\t";
      }
    }
    print '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg19&position=chr'.$rH_snp->{"chr"}.'%3A'.$rH_snp->{"position"}.'")';
    print "\t".$rH_snp->{"fullEffect"};
    print "\n";
  }
}

sub parseSampleNormalTumor {
  my $normalTumorFile = shift;
  my %normals;
  my %tumors;
  my %alias;

  open(SAMPLE, "$normalTumorFile") or die "Can't open sample normal tumor file: $normalTumorFile\n";
  while(my $line =<SAMPLE>) {
    chomp $line;
    if($line =~ /^#/) {
      next;
    }
    my @values = split(/,/, $line);
    $normals{$values[1]} = $values[0];
    $tumors{$values[2]} = $values[0];
    $alias{$values[0]} = [$values[1],$values[2]];
  }
  close(SAMPLE);
  print STDERR "STEP 3: Load paired csv file\n";
  return (\%normals, \%tumors, \%alias);
}

my $currentKGLine;
#sub init1kgMAF {
#  my $kgFh = shift;
#  $currentKGLine = <$kgFh>;
#  chomp($currentKGLine);
#  print STDERR "STEP 1: Read-in 1000G MAF file\n";
#}

#sub add1kgMAF {
#  my $kgFh = shift;
#  my $rHoH_snp = shift;
#  my $key = shift;
#  my $rH_snp = $rHoH_snp->{$key};
#  while(1) {
#    my @values = split(/\t/, $currentKGLine);
#    my $chromosome = $values[0];
#    my $position = $values[1];

#    if($chromosome ne $rH_snp->{"chr"}) {
#      $currentKGLine = <$kgFh>;
#      if(!defined($currentKGLine)) {
#        last;
#      }
#      chomp($currentKGLine);
#      next;
#    }

#    if($position < $rH_snp->{"position"}){
#      $currentKGLine = <$kgFh>;
#      if(!defined($currentKGLine)) {
#        last;
#      }
#      chomp($currentKGLine);
#      next;
#    }

#    if($position > $rH_snp->{"position"}){
#      last;
#    }
#    else {
#      $rH_snp->{"1kg_global_maf"} = $values[4];
#      $rH_snp->{"1kg_asian_maf"} = $values[5];
#      $rH_snp->{"1kg_american_maf"} = $values[6];
#      $rH_snp->{"1kg_african_maf"} = $values[7];
#      $rH_snp->{"1kg_european_maf"} = $values[8];
#      last;
#    }
#  }
#}
