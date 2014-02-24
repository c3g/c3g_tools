#!/usr/bin/perl

use Time::localtime;
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
    push(@chromosomes, $_);
  }
  close(CHR);
  $currentChr=shift(@chromosomes);
}

my %cohort;
my %cohort_count;

sub parseCohortFile {
  my $cohortFile = shift;
  
  open(COHORT, "<", $cohortFile) or die "Can't open $cohortFile\n";
  while ( <COHORT> ){
        chomp($_);
	my @split = split( "\t", $_);
        $cohort{$split[0]} = $split[1];
	$cohort_count{$split[1]}++;
  }
  close(COHORT);
}

sub initVcfFiles {
  my $rA_vcfFiles = shift;  
  my $rH_samples = shift;

  my %fileDetails;
  my %samples;
  for my $vcfFile (@{$rA_vcfFiles}) {
    my ($fileName) = $vcfFile =~ /(.*)\.vcf/;
    #print STDERR $fileName."\n";
    my $fh;
    open($fh, "<", $vcfFile) or die "Can't open $vcfFile\n";
    my $line;
    while( $line = <$fh> ) { 
      if($line =~ /^#CHR/) {
        last;
      }
    }
    chomp($line);
    my @values = split(/\t/,$line);
    
    if(@values < 10) {
      warn "File $vcfFile doesn't have enough columns\n";
      close($fh);
      next;
    }

    $line = <$fh>;
    chomp($line);
    $fileDetails{'file'} = {};
    $fileDetails{'file'}->{'fh'} = $fh;
    $fileDetails{'file'}->{'currentLine'} = $line;

    my $length = scalar @values;

    for (my $index = 9 ; $index < $length ; $index++ ){
	
    	my $sample = $rH_samples->{ $values[$index] };
	#print STDERR $index."\t".$values[$index]."\t".$sample."\n";

    	if(!defined($sample)) {
      		warn "VCF Sample name and sample sheet don't match: ".$sample."\n";
      		next;
    	}
       
    	$fileDetails{'sample'}->{$values[$index]} = $index;
    }
    
  }
  return \%fileDetails;
}

sub nextSamples {
  my $rHoH_samplesFiles = shift;

  my $line = $rHoH_samplesFiles->{'file'}->{'currentLine'};
  my @lineValues = split(/\t/, $line);
  
  my $chromosome = $lineValues[0];
  my $position = $lineValues[1];

  my $interval = $chromosome."-".$position;

  my %samplesToPrint;
  for( my $idx = 2 ; $idx < 9 ; $idx++ ) {
	push( @{ $samplesToPrint{'varInfo'}{$interval} }, $lineValues[$idx] );
	#print STDERR $lineValues[$idx]."\n";
  }

  for my $sample ( keys(%{ $rHoH_samplesFiles->{'sample'} } ) ) {        
      $samplesToPrint{'sample'}{$interval}{$sample} = $lineValues[$rHoH_samplesFiles->{'sample'}->{$sample}] ;
      #print STDERR $sample."\t".$lineValues[$rHoH_samplesFiles->{'sample'}->{$sample}]."\n";
  }

  my $fh = $rHoH_samplesFiles->{'file'}->{'fh'}; 
  my $line = <$fh>;

  if( !defined($line) ) {
    $samplesToPrint{'exit_status'} = 0;
    if( defined($fh) ) {
      close($fh);
    }
  }

  else{
	chomp($line);
	$samplesToPrint{'exit_status'} = 1;
  	$rHoH_samplesFiles->{'file'}->{'currentLine'} = $line;
  
  }
  return \%samplesToPrint;
}

sub parseSheet {
  my $fileName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";

# Skip header

  <SAMPLE_SHEET>;
  while(<SAMPLE_SHEET>) {
    chomp;
    my $line = $_;
    $line =~ s/"//g;
    my @values = split(/,/, $line);
    if($values[15] =~ /invalid/) {
      warn "Invalid: $values[0] $values[9] $values[11]\n";
      next;
    }

    my $runId = $values[9];
    #if($runId =~ /^0*([1-9]\d*)$/) {
    #  $runId = $1;
    #}

    my %sampleInfo;
    $sampleInfo{'name'} = $values[0];
    $sampleInfo{'libraryBarcode'} = $values[3];
    $sampleInfo{'runId'} = $runId;
    $sampleInfo{'lane'} = $values[11];
    #print STDERR "$values[0] $values[9] $values[11]\n";
    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

sub main {
  #my ($chromosomeFile, $kgFile, $sampleSheetFile, $cohortFile, @vcfFiles) = @ARGV;
  my ($chromosomeFile, $sampleSheetFile, @vcfFiles) = @ARGV;

  parseChr($chromosomeFile);

  my $rA_SampleLaneInfos = parseSheet($sampleSheetFile);

  my %sampleInfo;
  for my $rH_Sample (@$rA_SampleLaneInfos) {;
    if(!defined $sampleInfo{ $rH_Sample->{'name'} }) {
      $sampleInfo{ $rH_Sample->{'name'} } = [];
    }

    push(@{$sampleInfo{ $rH_Sample->{'name'} }}, $rH_Sample);
  }

  #parseCohortFile($cohortFile);

  my $rHoH_samplesFiles = initVcfFiles(\@vcfFiles, \%sampleInfo);
  #exit;
  my @samplesGiven = keys( %{ $rHoH_samplesFiles->{'sample'} } );
  #print STDERR join(" ",@samplesGiven)."\n";
 	print "chromosome\tposition\tref_allele\talt_alleles\tmappability_flag\tgene_name\tNb_Samples\tTotal_Depth\tFrac_RefDepth\tFrac_AltDepth\tAlt_QUAL\tAve_GenoQual\tMapQual\tDBSNP_GMAF\t";

        #foreach my $cohort ( sort { $a cmp $b } keys %cohort_count  ){
	#	print "$cohort\_count\tfrac_$cohort\t$cohort\_sample_id\t"

	#}

        #print "gene_description\tgo_ids\tgo_terms\t
	print "impact\teffect\teffect_type\tcodon_change\tamino_acid_change\tgene_type\ttranscript\tdb_snp_id\t";
	#print "1000genome_population_mafs\t";
 	print "CosmicID\tInterPro_Domain\tUniprot\tGERP++_neutral_rate\tGERP++_RS_score\t29way_SiPhy_score\tPolyphen2_prediction_based_on_HumVar_'D'_('probably_damaging')_'P'_('possibly_damaging')_'B'_('benign')\tSIFT_score\t";
 	print "Normal HomRef\tNormal HomAlt\tNormal Het\t";
 	for my $sampleName (@samplesGiven) {
   	  print "Depth_".$sampleName."\t";
   	  #print "GQ_".$sampleName."\t";
	  print $sampleName."\t";
 	}
 	print "ucsc_url\tfull_effect_description\n";


  while(1) {
    	my $rHoA_samplesToPrint = nextSamples($rHoH_samplesFiles);
    	if( $rHoA_samplesToPrint->{'exit_status'} == 0 ) {
      		exit;
    	} 
	else {
     
    		my $rHoH_snps;
    		for my $interval (keys(%{$rHoA_samplesToPrint->{'sample'}})) {
      			#$rHoH_snps = getSnpInformation($interval, $kgFh, $rHoA_samplesToPrint);
      	   		$rHoH_snps = getSnpInformation($interval, $rHoA_samplesToPrint);
    		}
    
  	#  printSnps(\@samplesGiven, \%snps, \%cohort);
    	printSnps(\@samplesGiven, $rHoH_snps);
	}
  }
}

my %sample_genotype;

sub getSnpInformation {
  my $interval = shift;
  my $rHoA_samplesToPrint = shift;
  #my $rHoH_snps = shift; 

  my @interval = split("-", $interval); 

  my $chromosome = $interval[0];
  my $position = $interval[1];

  my $rHoH_snps;
  my %snp;
  $snp{"chr"} = $chromosome;
  $snp{"position"} = $position;
  $snp{"dbSnpID"} = $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[0];
  $snp{"refAllele"} = $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[1];
  $snp{"altAllele"} = $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[2];
  $snp{"qual"} = $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[3];
   
  my $fullEffect;
  my $genotypeQuality;
  my $coverage;
  my $nmbRef;
  my $nmbAlt;
  my $MQ;
 
  my @readCountsForSnp;
  my %effects;

  $effects{"HIGH"} = [];
  $effects{"MODERATE"} = [];
  $effects{"LOW"} = [];
  $effects{"OTHER"} = [];

  my @infos = split(/;/, $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[5]);
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
    elsif($info =~ /^RO=(.*)/) {
      $snp{"RO"} = $1;
    } 
    elsif($info =~ /^AO=(.*)/) {
      $snp{"AO"} = $1;
    }
    elsif($info =~ /^DP=(.*)/) {
	$snp{"DP"} = $1;
      	#@readCountsForSnp = split(/,/,$1);
	#$snp{"DP4"} = \@readCountsForSnp;
    }
    elsif($info =~ /^DP4=(.*)/) {
      @readCountsForSnp = split(/,/,$1);
      $snp{"DP4"} = \@readCountsForSnp;
    }    

    elsif($info =~ /^MIL=(.*)/) {
      $snp{"mappabilityFlag"} = $1;
    }
    elsif($info =~ /^MQ=(.*)/){
	$snp{"MQ"} = $1;
    }
    elsif($info =~ /^GMAF=(.*)/) {
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
  }  ##Variant INFO column 

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
  $snp{"bioType"} = $rH_defaultEffect->{"bioType"};
  $snp{"transcript"} = $rH_defaultEffect->{"transcript"};
  $snp{"fullEffect"} = $fullEffect;

  #if($infos[0] =~ /INDEL/) {
  #  $interval .= "-INDEL";
  #}

  if(!defined($rHoH_snps->{$interval})) {
    $snp{"samples"} = {};
    $rHoH_snps->{$interval} = \%snp;
  }
  else {
    if(!($infos[0] =~ /INDEL/)) {
      if($rHoH_snps->{$interval}->{"refAllele"} ne $snp{"refAllele"}) {
        warn "Ref alleles don't match: $interval\n";
      }

      my @alleles = keys(%{ $snp{"altAllele"} });
      my $altAllele = $alleles[0];
      if(!defined( $rHoH_snps->{$interval}->{"altAllele"}->{"$altAllele"}) ) {
	#warn "Info: ".$infos[0]."\n";
	#warn "Alt alleles don't match: $interval, $sampleName\n";
        $rHoH_snps->{$interval}->{"altAllele"}->{"$altAllele"} = 1;
      }
    }
  }

    my @sampleFormat = split(/:/, $rHoA_samplesToPrint->{'varInfo'}->{$interval}->[6]);

    my $GQidx=-1;
    my $GTidx=-1;
    my $DPidx=-1;
    my $PLidx=-1;
    my $SPidx=-1;

    for(my $idx=0 ; $idx < @sampleFormat; $idx++) {
       if($sampleFormat[$idx] eq "GQ") {
         $GQidx=$idx;
       }
       if($sampleFormat[$idx] eq "GT") {
	 $GTidx=$idx;
       }
       if($sampleFormat[$idx] eq "DP") {
         $DPidx=$idx;
       }
       if($sampleFormat[$idx] eq "PL") {
         $PLidx=$idx;
       }
       if($sampleFormat[$idx] eq "SP") {
         $SPidx=$idx;
       }
    }

  for my $sample ( keys %{ $rHoA_samplesToPrint->{'sample'}->{$interval} } ) {

   	my @sampleInfo = split( /:/, $rHoA_samplesToPrint->{'sample'}->{$interval}->{$sample} );   

	$rHoH_snps->{$interval}->{$sample}->{"metrics"} = $rHoA_samplesToPrint->{'sample'}->{$interval}->{$sample};
	$rHoH_snps->{$interval}->{$sample}->{"DP"} = $sampleInfo[$DPidx];
  	$rHoH_snps->{$interval}->{$sample}->{"GQ"} = $sampleInfo[$GQidx];
	$rHoH_snps->{$interval}->{$sample}->{"GT"} = $sampleInfo[$GTidx];  
  	$rHoH_snps->{$interval}->{$sample}->{"PL"} = $sampleInfo[$PLidx];
	$rHoH_snps->{$interval}->{$sample}->{"SP"} = $sampleInfo[$SPidx];
  
	#print STDERR $sampleInfo[$DPidx]."\t".$sampleInfo[$GQidx]."\t".$sampleInfo[$GTidx]."\t".$sampleInfo[$PLidx]."\t".$sampleInfo[$SPidx]."\n";

  	if($sampleInfo[$GTidx] ne '0/0') {
    		if(substr($sampleInfo[$GTidx],0,1) ne substr($sampleInfo[$GTidx],2,1)) {
      			$rHoH_snps->{$interval}->{"samplesHet"} = $rHoH_snps->{$interval}->{"samplesHet"} + 1;
      			push( @{$sample_genotype{$interval}{"samplesHet"}},$sample );
    		} 
    		else {
      			$rHoH_snps->{$interval}->{"samplesHomAlt"} = $rHoH_snps->{$interval}->{"samplesHomAlt"} + 1;
      			push(@{$sample_genotype{$interval}{"samplesHomAlt"}},$sample);
    		}
  	}
  	else {
	 	$rHoH_snps->{$interval}->{"samplesHomRef"} = $rHoH_snps->{$interval}->{"samplesHomRef"} + 1;
	  	push(@{$sample_genotype{$interval}{"samplesHomRef"}},$sample);
  	}
  
  	#my @alleles = keys(%{ $snp{"altAllele"} });
  	#my @trueAlleles = split(/,/, $alleles[0]);

  	#my $newGenotype = $sampleInfo[$GTidx];
  	#$newGenotype =~ s/0/$snp{"refAllele"}/g;
  	#for(my $idx=0; $idx < @trueAlleles; $idx++) {
    	#	my $off = $idx+1;
    	#	$newGenotype =~ s/$off/$trueAlleles[$idx]/g;
  	#}
  	#$rHoH_snps->{$interval}->{$sample}->{"GT"} =~ s/$sampleInfo[$GTidx]/$newGenotype/g;
  }
  return( $rHoH_snps ); 
  #add1kgMAF($kgFh, $rHoH_snps, $key);
}

sub printSnps {
  my $rA_samplesGiven = shift;
  my $rHoH_snps = shift;
  #my $rH_cohort = shift;

  for my $rH_snp ( values( %$rHoH_snps ) ) {
    my $sampleSection = "";
    my $nbSamples     = 0;
    my $totalDepth    = 0;
    my $refDepth      = 0;
    my $altDepth      = 0; 

    my $GQ_sum = 0;

    my $key = $rH_snp->{"chr"}.'-'.$rH_snp->{"position"};

    for my $samples (sort @{$rA_samplesGiven}) {
	if( defined($rH_snp->{$samples}) ) {
		$GQ_sum += $rH_snp->{$samples}->{"GQ"};		
		#print STDERR $samples."\t".$rH_snp->{$samples}->{"GQ"}."\t".$GQ_sum."\n";

        }
    }

#    my %sample_cohort = ();

#     for my $samples (sort @{$rA_samplesGiven}) {
# 	if( defined($rH_snp->{"samples"}->{$samples}) ) {
# 	    if( defined( $rH_cohort->{$samples}  ) ){
# 		#print STDERR $samples."\n";
# 		push( @{ $sample_cohort{ $rH_cohort->{$samples} } }, $samples );	   
# 
#             }
# 
#         }
#     }

    if( defined ( $rH_snp->{"DP4"} ) ) {
          $totalDepth += $rH_snp->{"DP4"}->[0] + $rH_snp->{"DP4"}->[1] + $rH_snp->{"DP4"}->[2] + $rH_snp->{"DP4"}->[3];
          $refDepth   += $rH_snp->{"DP4"}->[0] + $rH_snp->{"DP4"}->[1];
          $altDepth   += $rH_snp->{"DP4"}->[2] + $rH_snp->{"DP4"}->[3];
    }

    if($totalDepth == 0){
    next;
    }
      
     
    print $rH_snp->{"chr"};
    print "\t".$rH_snp->{"position"};
    print "\t".$rH_snp->{"refAllele"};
    print "\t".$rH_snp->{"altAllele"};
    print "\t".$rH_snp->{"mappabilityFlag"};
    print "\t".$rH_snp->{"gene"};
    print "\t".@{$rA_samplesGiven};
    #print "\t\t\t";
    #printf ("\t%.1f",($nbSamples/@{$rA_samplesGiven}));
    print "\t".$totalDepth;
    #printf ("\t%.3f",($rH_snp->{"RO"}/$rH_snp->{"DP"}));
    #printf ("\t%.3f",($rH_snp->{"AO"}/$rH_snp->{"DP"}));
    printf ("\t%.3f",($refDepth/$totalDepth));
    printf ("\t%.3f",($altDepth/$totalDepth));

    print "\t".$rH_snp->{"qual"};
    printf ("\t%.1f",($GQ_sum/@{$rA_samplesGiven}));
    print "\t".$rH_snp->{"MQ"};
    if($rH_snp->{"GMAF"} eq "") {
                $rH_snp->{"GMAF"} = -1;
    }
    print "\t".$rH_snp->{"GMAF"};

#    foreach my $cohort (sort { $a <=> $b } keys %cohort_count){
#	if( defined( $sample_cohort{$cohort}  ) ) {
#    		print "\t".@{$sample_cohort{$cohort}};
#    		printf ("\t%.3f",(@{$sample_cohort{$cohort}}/$cohort_count{$cohort}));
#        	print "\t".join( ';',@{$sample_cohort{$cohort}} );
#	}else{
#		print "\t0\t0\t0";
#	}
#    }
#    print "\t".$rH_snp->{"geneDescription"};
#    print "\t".$rH_snp->{"goIds"};
#    print "\t".$rH_snp->{"goTerms"};
    print "\t".$rH_snp->{"impact"};
    print "\t".$rH_snp->{"effect"};
    print "\t".$rH_snp->{"effectType"};
    print "\t".$rH_snp->{"codonChange"};
    print "\t".$rH_snp->{"aaChange"};
    print "\t".$rH_snp->{"bioType"};
    print "\t".$rH_snp->{"transcript"};
    print "\t".$rH_snp->{"dbSnpID"};
#    print "\tASN:".$rH_snp->{"1kg_asian_maf"}.' AMR:'.$rH_snp->{"1kg_american_maf"}.' AFR:'.$rH_snp->{"1kg_african_maf"}.' EUR:'.$rH_snp->{"1kg_european_maf"};
    print "\t".$rH_snp->{"cosmic"};
    print "\t".$rH_snp->{"interPro"};
    print "\t".$rH_snp->{"uniprot"};
    print "\t".$rH_snp->{"gerpNR"};
    print "\t".$rH_snp->{"gerpRS"};
    print "\t".$rH_snp->{"29way"};
    print "\t".$rH_snp->{"polyphenHvar"};
    print "\t".$rH_snp->{"sift_score"};
    print "\t".(0+$rH_snp->{"samplesHomRef"});

    if(defined( $sample_genotype{$key}{"samplesHomRef"} ) ){
    	    print "(".join( ';',@{$sample_genotype{$key}{"samplesHomRef"}} ).")";
    }
    print "\t".(0+$rH_snp->{"samplesHomAlt"});
    if(defined( $sample_genotype{$key}{"samplesHomAlt"} ) ){
	    print "(".join( ';',@{$sample_genotype{$key}{"samplesHomAlt"}} ).")";
    }
    print "\t".(0+$rH_snp->{"samplesHet"});
    if(defined( $sample_genotype{$key}{"samplesHet"} ) ){
    	    print "(".join( ';',@{$sample_genotype{$key}{"samplesHet"}} ).")";
    }
    %sample_genotype = ();
    print "\t";
   	for my $sampleName (@{$rA_samplesGiven}) {
   	  if( defined( $rH_snp->{$sampleName} ) ) {
   	    print $rH_snp->{$sampleName}->{"DP"}."\t";
	    #print $rH_snp->{$sampleName}->{"GQ"}."\t";
   	    print $rH_snp->{$sampleName}->{"metrics"}."\t";
   	  }
      else {
        print "\t";
        print "\t";
      }
    }
    print '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg19&position=chr'.$rH_snp->{"chr"}.'%3A'.$rH_snp->{"position"}.'")';
    print "\t".$rH_snp->{"fullEffect"};
    print "\n";
  
  }
}

my $currentKGLine;
sub init1kgMAF {
  my $kgFh = shift;
  $currentKGLine = <$kgFh>;
  chomp($currentKGLine);
}

sub add1kgMAF {
  my $kgFh = shift;
  my $rHoH_snp = shift;
  my $key = shift;
  my $rH_snp = $rHoH_snp->{$key};
  while(1) {
    my @values = split(/\t/, $currentKGLine);
    my $chromosome = "chr".$values[0];
    my $position = $values[1];

    if($chromosome ne $rH_snp->{"chr"}) {
      $currentKGLine = <$kgFh>;
      if(!defined($currentKGLine)) {
        last;
      }
      chomp($currentKGLine);
      next;
    }

    if($position < $rH_snp->{"position"}){
      $currentKGLine = <$kgFh>;
      if(!defined($currentKGLine)) {
        last;
      }
      chomp($currentKGLine);
      next;
    }

    if($position > $rH_snp->{"position"}){
      last;
    }
    else {
      $rH_snp->{"1kg_global_maf"} = $values[4];
      $rH_snp->{"1kg_asian_maf"} = $values[5];
      $rH_snp->{"1kg_american_maf"} = $values[6];
      $rH_snp->{"1kg_african_maf"} = $values[7];
      $rH_snp->{"1kg_european_maf"} = $values[8];
      last;
    }
  }
}
