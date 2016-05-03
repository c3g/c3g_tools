28 tags, 274 commits

HEAD        Tue Jan 26 10:55:10 2016 -0500        8 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       354ce82 corrected a bug in chipSeqgenerateAnnotationGraphs.R that caused a blank exon annotation graph - BFXDEV-501
       b0c47cc mugqic_tools - up-to-date CHANGELOG.md
       ab85318 ading a changelog to mugqic_tools (since v2.1.5)

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       96eadd0 R-tools - update SNV metrics to the new stats output fromat from snpEff - BFXDEV-490
       99b25f3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       e45cf73 Rtools - remove typo and allow to skip insert size when working on single end library - BFXDEV-502
       b7a6767 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       6863562 perltools - remove # character from blast hit results - BFXDEV-500

2.1.4        Fri Jan 8 12:28:52 2016 -0500        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       7d49d89 adding of ampliconseq tool

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       4c1e35d added script to generate gemini compatible vcfs

2.1.3        Wed Nov 25 16:43:54 2015 -0500        3 commits

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       12e7e61 Removed SNVstats error with missing lat line in the list file - BFXDEV-120
       c0f12ba Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       653c15d Adding python-tools/mergePairedSingleTrimFastq.py a tools to interlace paired fastqs and also concatenate single fastq in one majopr fastq after the trimming

2.1.2        Tue Sep 22 11:32:36 2015 -0400        16 commits

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      1 commits

       84ff003 pacBioAssemblyStats.pl edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       f2710d9 Fixed bug where annotations could not contain # char

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      11 commits

       bbdd52b BFXDEV-397 filtered tab separated assembly file extension must be tsv not xls
       bb8b1ee BFXDEV-397 join filtered isoforms and isoform counts matrix to generate exploratory analysis on filtered transcripts
       2ce156d BFXDEV-397 filter components using trinotate annotations generate a fasta and tab separated filtered assembly file
       2042d52 BFXDEV-397 extract sequences from assembly using python SeqIO, generate tabular and fasta filtered assembly
       a6d8d13 BFXDEV-397 added functionnality to filter files using python expressions, applied to parse Trinotate Output and to parse and merge csv files
       6b44025 BFXDEV-396 filter trinotate output to create a file with IDs of annotated transcripts given a set of rules in the ini file, added to parseMergeCsv the functionnality to rename columns that have the same name in input files (R make.names like)
       444b013 BFXDEV-396 parse merge csv, corrected errors in sort by when join gets None
       65337cf BFXDEV-396 parse trinotate output, read isoforms length and print only the longest isoform per gene in annotations, added a standard parse plus merge csv files python program
       08d2c29 PRJBFX-1075 mugqic_pipeline tests, scripts used to setup automatically a test directory (sync data and programs)
       4b82c77 BFXDEV-396 parse trinotate output to obtain gene annotations, transcript annotations and goterms
       65055cb corrected issue in validateBAMContent.pl when any field in the sample sheet contains a comma BFXDEV-406

  lletourn <louis.letourneau@mail.mcgill.ca>      2 commits

       187363f BFXDEV-367 added year directory flexibility
       52d8592 BFXDEV-370 Added a way to generate split bed files from a dictionary file

  ptranvan <patrick.tranvan@mail.mcgill.ca>      1 commits

       c6d9982 BFXDEV-434 Script for amplicon-seq pipeline.

2.1.1        Wed Mar 25 12:14:14 2015 -0400        3 commits

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       64b89fa changed threshold test to allow reporting FDR of 1

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       5d67c83 Fixed bug, array was interpreted instead of using @ as a character

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       fd4b011 python - rRNacounter.py - correct issues when missing transcript Id or Name in the gtf and add a shorten version of the stat file <OUTPUT>_short.tsv

2.1.0        Mon Mar 9 15:58:31 2015 -0400        5 commits

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       b8cd478 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       556d49d BFXDEV-59 DNAsampleMetrics.R fix

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       1fe0c67 python - update rrnaBAMcounter.py to include all unmapped reads in the total count

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       95a1618 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       721a2b3 python - add rrnaBAMcounter.py a small tools to report rrna counts

2.0.3        Fri Jan 30 10:55:17 2015 -0500        9 commits

  Joël Fillon <joel.fillon@mcgill.ca>      1 commits

       3bce44e BFXDEV-292 Updated MACS2 file paths in R-tools/chipSeqgenerateAnnotationGraphs.R

  lletourn <louis.letourneau@mail.mcgill.ca>      4 commits

       79092c7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       8896743 Added output file name instead of redirection to allow stream buffering
       5266243 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       70fd054 File needs to exist in paired file to add columns

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       b38b144 TOOLS - update tools/gtf2tmpMatrix.awk to support commetn lines in the gtf

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       b11c11c Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       6d02639 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b115892 Tools -add runLumpy.sh using v0.2.9 of Lumpy and the python script for the mean and stdev

2.0.2        Tue Dec 16 16:43:51 2014 -0500        7 commits

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       3d45241 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       133314b Fixed line names of PacBio read count statistics in pacBioAssemblyStats.pl

  lletourn <louis.letourneau@mail.mcgill.ca>      5 commits

       bee0f73 Merged
       ec283d9 Added check for remaining bams
       b353c17 changed shebang
       394e526 Removed smtps code
       1c25d6a Added name test

2.0.1        Wed Dec 10 13:55:58 2014 -0500        1 commits

  Joël Fillon <joel.fillon@mcgill.ca>      1 commits

       3d37484 Fixed bug rpkmSaturation.R read.table with options quote='', comment.char=''

2.0.0        Fri Dec 5 16:50:16 2014 -0500        5 commits

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       e72d9d7 README.md edited online with Bitbucket
       c559d56 Moved getLogReport.pl nanuq2mugqic_pipeline.py into mugqic_pipelines/utils/

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       2798ba3 BFXDEV-285 Added BED download to nanuq2mugqic

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       e18c5e4 R-tools - goseq.R : remove conflict
       ac10aed R-tools - goseq.R : allows more flexibity in fdr and p-value for goseq AND remove native goseq approach - BFXDEV-289

1.10.5        Wed Nov 19 17:45:08 2014 -0500        5 commits

  Joël Fillon <joel.fillon@mcgill.ca>      3 commits

       df16cb8 Fixed read.table parsing with more options in goseq.R
       3cc7d88 Renamed readset and design file 'SampleID' column into 'Sample'
       66a0ef5 Minor cosmetic change in nanuq2mugqic_pipeline.py

  lletourn <louis.letourneau@mail.mcgill.ca>      2 commits

       0fe9969 fixed multi region sample problem.
       7f9d6f1 fixed multi region sample problem.

1.10.4        Mon Oct 27 10:14:56 2014 -0400        8 commits

  Joël Fillon <joel.fillon@mcgill.ca>      4 commits

       4bd19a3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       4c856bb Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       2895b24 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       69045a5 Change column name Sample -> SampleID

  lletourn <louis.letourneau@mail.mcgill.ca>      3 commits

       5b7edf2 Removed dependency
       248b569 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       0cabf7b Added configurable depth

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       6a5a015 gtf2tmpMatrix.awk: add a space to underscore replacement - BFXDEV-275

1.10.3        Fri Sep 19 14:26:49 2014 -0400        2 commits

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       f9e7384 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       3549781 R-tools/DNAsampleMetrics.R: allow dots in sample names - BFXDEV-265

1.10.2        Thu Sep 18 14:20:04 2014 -0400        21 commits

  eric audemard <eaudemar@imac6-ub.(none)>      2 commits

       5b255ea Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_tools
       8142961 fixe bug for puure

  Joël Fillon <joel.fillon@mcgill.ca>      5 commits

       b41f22a Fixed some wrong shebangs with #!/usr/bin/env ... + moved fasta_Length_select.py in python-tools/
       ff00af7 Added Pacbio support in python-tools/nanuq2mugqic_pipeline.py
       e20a00c Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       0b7d384 Minor nanuq path update
       79ea08e Added project ID in filenames created by nanuq2mugqic_pipeline.py

  lletourn <louis.letourneau@mail.mcgill.ca>      14 commits

       ef5b466 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b2df75a BFXDEV-257 fixed 'track' bed header handling
       102d25f Fixed sampleName source from validateBAMContent.pl
       9dcdfeb Added sampleSetup filters
       228da50 Fixed STARTTLS
       878bcac Fixed bad line
       9a9d8f3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       12514d2 Add smtp server setting to tool
       91925ce Fixed relative path bug
       dabdf6a Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       f5f2705 Added a tool to track BAM contents vs a sampleSheet
       e286095 Added a usage computation script
       f072cf0 BFXDEV-239 fetch bed files automatically
       d2f5483 Fixed bug

1.9.3        Wed Aug 13 15:40:58 2014 -0400        0 commits

1.10.1        Wed Aug 13 15:42:47 2014 -0400        23 commits

  eric audemard <eaudemar@imac6-ub.(none)>      3 commits

       50a26f2 add awk script for puure pipeline
       aec0ad3 rm tempory file
       f086553 add python and R script for PUURe

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       e47fc25 check.names=FALSE in rpkmSaturation.R to avoid sample names renaming

  Joël Fillon <joel.fillon@mcgill.ca>      5 commits

       10b2015 Minor help formatting in nanuq2mugqic_pipeline.py
       e04b1d1 Renamed nanuq2readset.py into nanuq2mugqic_pipeline.py
       0651d43 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools into python
       a3c469b Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools into python
       4e195ec Added nanuq2readset.py

  jtrembla <jtremblay514@gmail.com>      6 commits

       75b61db Added bam to fastq script.
       149c76d Added argument for low abundance cluster cutoff after 99% clustering step. BFXDEV-31
       0db65b3 Fixed bug when selecting Fungi in OTU tables for 18S analyses. BFXDEV-31
       91a3cb0 Minor fix for cluster counts in report table. BFXDEV-31
       a3977a6 Added 454 option to sampleSetup.pl. BFXDEV-31
       49576a9 Added script to keep best n blast hit from blast tables. BFXDEV-30

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      4 commits

       5f77252 fixed verbose level. BFXDEV-31
       bef496e Added --sampleSheet argument to getMiSeqBarcodes.pl . BFXDEV-31
       5ecef7e Improved script to merge 454 fastqs prior to do 16S/ITS analysis. BFXDEV-31
       705ebbe Refined support for 454 libs. Works with amplicon data at least... BFXDEV-31

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      1 commits

       aae7d76 Added naive interface script to gqSeqUtils::calculateGenelengthsFromGtf

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       03494cb perl-tools/csvDepthFilter.pl
       c2acd2c README.md edited online with Bitbucket
       9dea6af README.md edited online with Bitbucket

1.10        Thu Jun 5 15:20:49 2014 -0400        7 commits

  jtrembla <jtremblay514@gmail.com>      1 commits

       7c2b7cc BlastCov.tsv is now sorted. BFXDEV-30

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       bb3351c Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       5855f97 Add tools/WG2ChromosomeFasta.awk to split wholeGenome fatsa into chromosme fasta

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       f583128 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       66cf8b6 DNAsampleMetrics.R will add NA filed if missing file of CCDS metrics - BFXDEV-214
       7fc08bd R-tools/deseq.R and  R-tools/edger.R now support numeric sample name starting with 0 in the design file - BFXDEV-213

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      1 commits

       f18b168 genome length R script

1.9        Thu May 22 15:52:21 2014 -0400        32 commits

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       8dedf91 Reformatted FASTQ/BAM check in sampleSetup.pl
       36496aa sampleSetup with BAM/FASTQ support

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      1 commits

       90ca70d Automatic install genomes assumes that variable MUGQIC_INSTALL_HOME is already set, commented the prologue section

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      2 commits

       0bbfb7a Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       e761d84 BFXDEV-133 Bug detected in guillimin phase2 : LoadConfig must not be used outside  pipeline wrappers

  jtrembla <jtremblay514@gmail.com>      17 commits

       2e10fea Fixed indentations RRNATagger scripts. BFXDEV-31
       96a2cc3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       8c234b8 Modified percent cutoff to X cov. BFXDEV-30
       72beb43 Script to add sample names to a vcf file. Have to provide a sample list in a separate txt file.
       0c40368 Fixed XML file output. BFXDEV-30
       6427884 Added option to output the runId. BFXDEV-30
       05c9794 Added estimated coverage in one of the generated tables. BFXDEV-30
       e1b188b Added script to add coverage in blast table (pacbio pipeline). BFXDEV-30
       ca16e2b Addition of script to merge FPKMs from an RNASeq project.
       d1544d0 Wrote script to merge many raw count matrix (form the RNASeq pipeline). PRJBFX-699
       ad33ba8 Script to compile Ray results.
       73e465d R script to plot histogram of contig length (use with ray output).
       bda1da9 Minor fixes. BFXDEV-31
       67647b1 Cut pacbio reads to 2000 nt max. Because of crash in downstream qscore fastx step. BFXDEV-31
       83520b1 Added verbosity for debugging. BFXDEV-31
       623fa15 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       fdefc43 added script template to perl tools.

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      2 commits

       9d508fb Added script to preprocess pacbio reads to make them compatible with RRNATagger pipeline. BFXDEV-31
       ea38231 Modified pacBio samplesetup so it can output ccs reads as well. BFXDEV-31

  lletourn <louis.letourneau@mail.mcgill.ca>      8 commits

       f5bbfb7 BFXDEV-196 Added onTarget calculations
       2228327 Deprecated
       4d5c477 Added tool to add basefreq to vcf2csv generated file
       09a3ce7 Added CAF support
       86af2d2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       7a25aa5 Fixed sample data validation
       609e706 Added paired version of vcf2csv
       ce20505 Added check on columns for FASTQs and BAM files

1.8        Fri Mar 21 10:50:47 2014 -0400        24 commits

  Joël Fillon <joel.fillon@mcgill.ca>      6 commits

       60d4c2b Fixed minor mug in mergeTrimmomaticStat.R: don't raise error if metrics directory already exists
       272534a Fixed exit code parsing with string instead of digits, in case of negative value
       9deeadd Fixed exit code parsing on guillimin phase2 + added job status
       a7232fe Replaced shebang #!/usr/bin/perl with #!/usr/bin/env perl in all Perl scripts
       38e99ca Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       596ae25 Added bam support in sampleSetup.pl

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      2 commits

       f82a5ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       90e9ade corrected bugs found installing genomes on Guillimin

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      2 commits

       c847345 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       65393b6 Convert vcf generated by the standard pipeline (samtools+snpeff+snpsift annotations) to a tab separated file

  jtrembla <jtremblay514@gmail.com>      5 commits

       a3ffbfc Modified summary table in pacBioAssemblyStats.pl. BFXDEV-30
       6250902 Updates to pacbio scripts. BFXDEV-31.
       69d9699 Fixed font size. BFXDEV-31
       77f3a7e Added OTUheatmap.R for RRNATagger pipeline. BFXDEV-31
       5ffccd1 added plotTree.R for RRNATAgger pipeline. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      2 commits

       f6f56c0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       1938c7b tool to plot coverage bed files. BFXDEV-30

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       64d2d44 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       4edc18e R-tools/DNAsamplemetrics.R: change the alignment and duplicates metrics which where not good when mem is used correct BFXDEV-130
       8f31f92 remove unneccessary file integrity in vcfstats.py

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       fe831ce Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       90c9d92 chipSeqgenerateAnnotationGraphs.R : add additional control for empty file with size > 0
       3c5d163 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       6f4cab2 add null buffer in vcfstat.py; correct BFXDEV-120

1.7        Thu Feb 27 14:21:48 2014 -0500        22 commits

  Joël Fillon <joel.fillon@mcgill.ca>      1 commits

       11a876c Added count of successful/active/inactive/failed jobs in getLogReport.pl

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      6 commits

       adf0ba6 BFXDEV-118 corrected typo : Error in ccdswgbase500[nameLoc] = statsValue[[1]][12] : object 'ccdswgbase500' not found
       4ffbc51 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       8979f4a BFXDEV-82
       1e5bb9a BFXDEV-90 additional parameters could be added when downloading fasta sources. Use exonerate to create the byChr fasta files
       869df6b set and unset IFS to add extra parameters for wget
       6d2e5f2 add the byChr directory using the exonerate's fastaexplode tool

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      8 commits

       060424c Added sampling of 5000 reads instead of 1000 to determine if phred 64 or 33. BFXDEV-31
       ac0e6f6 Added tools associated for RRNATagger pipeline. BFXDEV-31
       1d195ad Script showing which perl modules are currently installed (based on your PERL5LIB env.).
       6103763 Tool to fetch MiSeq barcodes having the runId in hand. Generates barcodes fasta file for that run. BFXDEV-31
       0b13647 Added debug info in STDERR so we can "more easily" associate Xpercent values with X cov. values. BFXDEV-30
       d21120c Fixes for contigs N25, N50, etc. Support for multiple polishing rounds in stats ouput. BFXDEV-30
       a41ff66 Modified path in which script looks for basx files (i.e. /2014/..) BFXDEV-30
       877029f Updated scripts to compute stats from polished assembly. + added HGAP cutoff value in X (cov). BFXDEV-30

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       7137697 add a checking step for a correct add of the change rate file name at the end of vcfStats.py. because sometime it doesn't works and don't give exit code 1 - BFXDEV-120
       355ef9d Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       bcbc074 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       18dbf68 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       9dfc311 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       630ede6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b1da309 update R-tools/deseq.R round count to allow RSEM estimates

1.6        Thu Jan 16 11:50:26 2014 -0500        20 commits

  Joël Fillon <joel.fillon@mcgill.ca>      4 commits

       81edd0d Added option to keep unsuccessful jobs only
       31c1a92 Use array instead of hash for memtime values
       b0b7273 Added --success option + updated memtime lines with MEMTIME prefix
       1c1ca04 Added in log report lowest/highest memory jobs, queue, limits and many other columns

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      1 commits

       9b03a64 bug: bwa and bowtiw indexes in their respective directories

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      7 commits

       6301360 updated pacbio R script to generate plots. BFXDEV-30
       f8f88f5 Added one missing \t in header line when --memtime option. BFXDEV-30
       5284a79 Removed STDERR debug lines . BFXDEV-30
       fdf618d Creates symlinks in raw_reads now.
       1c99b73 Added support for memtime output in job output files. BFXDEV-30
       0b47d59 Added script that extracts and compile assembly stats from multiple assembly done with the pacbio assembly pipline.
       d3561c8 Updates to pacbio scripts for pacbiopipeline. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      2 commits

       40eb4ad BFXDEV-81 Force encoding of Text::CSV
       d6175d8 BFXDEV-60 added error handling in blast script

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       20bee87 remove integer formating edger
       be18ef7 remove iunteger formating deseq
       46cdbfe rounding and formating raw_count-matrix object to allow RSEM output
       92150f5 edger.R edited online with Bitbucket
       83a3709 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       2cff02b update R-tools/deseq.R R-tools/edger.R tools/gtf2geneSize.awk

1.5        Wed Dec 4 17:31:54 2013 -0500        9 commits

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       fae748f removed method='deviance' and robust=TRUE. Dispersion estimates otherwise crashed with large number of tested features.

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       36107a7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       97dd764 Added missing executable permissions on some scripts

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       b3c7b21 Remove die if missing xml files. BFXDEV-30
       3734d18 Fixed some display aspect of plots. BFXDEV-30
       5d38409 Plot data on filtered reads only. BFXDEV-30
       e2c59ba Now accepting bax.h5 files. BFXDEV-30
       21e25c1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       9afd05f Added coverageFraction options. BFXDEV-30

1.4        Tue Nov 26 16:15:09 2013 -0500        14 commits

  Joël Fillon <joel.fillon@mcgill.ca>      4 commits

       32bd762 Added parallelBlast.pl from David Morais
       8aa30c7 Added extra virtual memory percentage
       1f5c951 Minor fix
       e75b39c Improved text log report with human readable durations and memory size + ratios

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      10 commits

       baa2e84 Added parsing feature for xml file. BFXDEV-30
       124242b Removed print debug lines.. BFXDEV-30
       30df1a6 Added fastq support for hgap cutoff. BFXDEV-30
       d205d97 Fix because new columns where added on nanuq samplesheets. BFXDEV-30
       e55569b Fixes + compute reads stats directly form .tsv. BFXDEV-30
       1cc0b4b Skip $filteredSummaryTable if file does not exists. BFXDEV-30
       bcfd86d minor fix. BFXDEV-30
       f30195f Minor fix to plot script. BFXDEV-30
       af895cd Fixes to PacBio scripts. BFXDEV-30
       3154bf0 Skipped header line for DNAsampleMetrics.R. Plus a bunch of fixes for pacBio perl scrips...

1.3        Tue Nov 19 11:26:09 2013 -0500        11 commits

  Joël Fillon <joel.fillon@mcgill.ca>      1 commits

       ba0c093 getLogReport now handle job output relative paths

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      5 commits

       8d1ab10 chipseq annotations plots : changed graphs title-font size
       ab07c78 Print header on read/alignment metrics
       16ce33c chipseq annotations plots : changed graphs title
       d628bb6 Plots of annotation stats chipseq pipeline
       72db555 added parser for chipseq annotations - narrow peaks

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      2 commits

       e5924a1 Added sequence iterators.
       86b53d8 Added R script to generate pacbio reads distribution plots. More to come...

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       83fa74f Added a more verbose explanation of sampleSetup

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       83f5bed Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       885d006 correct vcfStat.py exception -  BFXDEV-50

1.2        Mon Nov 4 16:03:51 2013 -0500        11 commits

  Joël Fillon <joel.fillon@mcgill.ca>      2 commits

       59ae2f0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       4ba58a2 Moved former Log lib into getLogReport.pl

  jtremblay@genomequebec.com <jtrembla@abacus2.(none)>      1 commits

       e92437f Added EOL ($) for each patterns to not consider .mugqic.done files.

  lletourn <louis.letourneau@mail.mcgill.ca>      2 commits

       853f9ad Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       fb1406b Added library barcode support

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       4a29ecd Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       5bda081 update edger.R
       6a456ff R-tools/snvGraphMetrics.R update
       0ee0703 remove the double $ patern in DNAsampole Metrics
       ef7d3e1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       445a476 snvGraphicMetrics.R add robustness to single sample analysis

1.1        Thu Oct 24 14:02:21 2013 -0400        1 commits

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       c87cd72 correct typo in DNA Metrics script

1.0        Thu Oct 24 10:50:35 2013 -0400        1 commits

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       52f686f add regex to avoid .done captation and rename the tag from v1.0 to 1.

v1.0        Tue Oct 22 16:20:47 2013 -0400        3 commits

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      3 commits

       1902529 Rtools snvGraphMatreics.R for the dnaSeq pipeline
       7dfaa57 update snvGraphs
       bde9560 update SNV metrics

