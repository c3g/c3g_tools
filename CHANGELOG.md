81 tags, 908 commits

HEAD        Mon Jul 7 10:33:11 2025 -0400        4 commits

  MareikeJaniak <30602072+MareikeJaniak@users.noreply.github.com>      1 commits

       dd3ac91 Merge pull request #2 from c3g/rp_metrics_fixes

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      3 commits

       c298d35 Version bump to 2.12.12
       e9e9254 bug fix to choose between parsing of bcl2fast or bclconvert
       d8f85ed Version bump to 2.12.11

2.12.11        Fri Jul 4 11:06:51 2025 -0400        5 commits

  MareikeJaniak <30602072+MareikeJaniak@users.noreply.github.com>      1 commits

       a409baa Merge pull request #1 from c3g/rp_metrics_fixes

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      4 commits

       2da5082 fix zero coverage error for chrY
       478c312 remove comment
       902aacf fixes for bclconvert metrics, insert size metrics
       a1deb03 Version bump to 2.12.10

2.12.10        Fri Jun 6 09:13:41 2025 -0400        18 commits

  Jean-Michel Garant <jean-michel.garant@computationalgenomics.ca>      5 commits

       a0cc255 Fix list concatenation of fastp quality curves mean
       8262dff Merged in run_processing (pull request #41)
       35b370b Support single end sequencing metricstojson for Freezeman
       bab2c25 Nonetype instead of N/A string to comply with Freezeman ingestion
       61c35a9 Include median insert size to the run_validation_json

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      6 commits

       6f74e0a Version bump to 2.12.10
       35a1715 update release instructions for github
       e708592 Merge branch 'master' of github.com:c3g/c3g_tools
       5449cdd formatting
       c0812c7 updating run processing metrics parsers for bclconvert
       3fac7e7 in run processing set sex_match to None if sex_det is equal to ? to avoid false sex concordance FALSE in multiqc report

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      6 commits

       a05f8d4 Merged in bcl_convert_parsing (pull request #42)
       893ddf9 Merged in rp_blast (pull request #40)
       01e4a40 limit parallel jobs for runBlast.sh
       39e1e14 replace standard blast with parallel blast for run processing
       d232016 Version bump to 2.12.10-beta
       57e6967 Version bump to 2.12.9

  Paul STRETENOWICH <31796146+paulstretenowich@users.noreply.github.com>      1 commits

       1b9ac43 Update README.md

2.12.9        Thu Oct 17 11:11:42 2024 -0400        5 commits

  Jean-Michel Garant <jean-michel.garant@computationalgenomics.ca>      3 commits

       c769320 Merged in run_processing (pull request #38)
       ced1b26 Fix the math for inferred sex
       e645189 Sex concordance value is can be set to None/null

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      2 commits

       3d731db Version bump to 2.12.9-beta
       015f3cb Version bump to 2.12.8

2.12.8        Thu Jul 11 11:13:14 2024 -0400        4 commits

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      2 commits

       7a86718 Version bump to 2.12.8-beta
       453fee4 Version bump to 2.12.7

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      2 commits

       42c12e9 Merged in tools_for_GenPipesV (pull request #37)
       48863c6 add script to prepare kallisto report for multiqc

2.12.7        Wed May 15 13:49:40 2024 -0400        4 commits

  Jean-Michel Garant <jean-michel.garant@computationalgenomics.ca>      2 commits

       bcb9856 Merged in metrics_types (pull request #36)
       8a7d6b3 Add checks to make sure alignment metrics are passed as numerals not strings

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      2 commits

       b3146b5 Version bump to 2.12.7-beta
       34d6ab7 Version bump to 2.12.6

2.12.6        Fri Dec 22 14:01:37 2023 -0500        5 commits

  Jean-Michel Garant <jean-michel.garant@computationalgenomics.ca>      1 commits

       39be23f Removing stray commas that causes tuples & arrays in json

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      3 commits

       7907307 python_tools : fastp metrics to rp json
       98eb063 Version bump to 2.12.6-beta
       2ec2d71 Version bump to 2.12.5

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      1 commits

       ccaa719 Merged in fastp_parser (pull request #35)

2.12.5        Wed Nov 29 10:09:37 2023 -0500        5 commits

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      2 commits

       93da164 python-tools : permissions
       5ed151d python tools : add splitBarcode parser

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      3 commits

       2c13ebe Merged in splitbarcode_parse (pull request #34)
       142f566 Version bump to 2.12.5-beta
       e1cbb84 Version bump to 2.12.4

2.12.4        Mon Oct 23 09:37:41 2023 -0400        9 commits

  Eric Audemard <eric.audemard@gmail.com>      2 commits

       b800e80 Merged in master (pull request #32)
       6ec622a Write output before generate figure. Write all peaks results in the output file (not only significant one)

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      3 commits

       e0d55c6 DiffBind : keep th adjustable parameter
       e3202a9 Version bump to 2.12.4-beta
       c2573df Version bump to 2.12.3

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      4 commits

       b4b5f9e Merged in diffBind_th (pull request #33)
       2be549b DiffBind : display plot errors in output instead of exiting
       780afc0 DiffBind : add th parameter to fig 3.2
       f764c87 diffBind : add second th parameter

2.12.3        Wed Oct 18 14:36:25 2023 -0400        6 commits

  Henrion <edouard.henrion@computationalgenomics.ca>      1 commits

       4924acd Version bump to 2.12.2

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      4 commits

       5d71c5a Version bump to 2.12.3
       42834ec mugqic_tools : merge changes from rjme
       890f006 run processing : add file sizes for fastqs and bams to run validation report
       093364a run processing : fixed error in check for input file exist

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      1 commits

       bc763b0 Merged in run_proc_file_sizes (pull request #31)

2.12.2        Thu Jun 8 21:24:45 2023 -0400        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       504679d Python Tools - RunProcessing : corrected bcl2fastq metrics

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       7720091 Version bump to 2.12.1

2.12.1        Mon May 29 14:33:36 2023 +0000        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       7988963 Version bump to 2.12.0

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       4956ae8 R-tools - fix typo in hicrep.R

2.12.0        Tue May 23 13:27:55 2023 -0400        7 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       0917db9 Merge remote-tracking branch 'remotes/origin/master' into HEAD

  ehenrion <edouard.henrion@mcgill.ca>      4 commits

       495d14e Merged in fix_asva_plot_error (pull request #30)
       1836eda Merged run_proc_match_undetermined into master
       5898e5d R Scripts - corrected hicrep.R for "no down-sampling" cases
       0557e81 Version bump to 2.11.3

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      1 commits

       43338f2 python tools : run_processing - match undetermined barcodes

  P-O Quirion <pioliqui@gmail.com>      1 commits

       6b62ccc bugfix cleanOTUtable

2.11.3        Fri May 5 17:52:27 2023 +0000        6 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       2df3ffd Python Tools - Updated getFastaBinedGC for Python3 compatibility
       cc9910b Version bump to 2.11.2

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      4 commits

       840ad06 Merged in diffBind_bam (pull request #29)
       d626d7c R-tools DiffBind : changed default bam file name
       29b6267 R-tools DiffBind : fixed typo
       e6afe40 R-tools DiffBind : added input option to specify extension of bam files

2.11.2        Tue Apr 4 12:53:12 2023 -0400        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       d0392e8 Improvements in runprocssing tools
       374d0a3 Version bump to 2.11.1

2.11.1        Fri Dec 2 10:47:35 2022 -0500        4 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       47765a1 Python tools - python-tools/AmpliconSeq_script.py : minor typo fix

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      1 commits

       4c7e2d7 CHANGELOG for 2.11.0

  ehenrion <edouard.henrion@mcgill.ca>      2 commits

       568c9cd FastqReadNameEdit updated
       ee3dc4a Some minor changes to R-tools

2.11.0        Wed Nov 9 15:16:49 2022 +0000        19 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      6 commits

       8278f1c undo nundesired changes in ihec_metrics
       ef03f1e minor format change
       fefa6a0 Python Scripts - updating run Processing metrics
       5d3ec55 Python Scripts - added working_dir as 2nd argument in ihec_metrics
       69da92f R-Scripts - Added batch efet corection to deseq2
       3f10c69 Version bump to 2.10.10

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      13 commits

       57b523a Version bump to 2.11.0
       920057b R-tools - correct typo in deseq2.R
       8913947 R-tools - corretions in deseq2.R
       be3b365 R-TOOLS - cleaning deseq2.R from debug prints
       da205b8 R-tools - some adjustments to the batch effect correction in deseq2.R
       0c180c5 minor
       fb9edf6 R-tools - Fixed deseq2.R for batch effect correction
       02599f9 Python-tools - updates in AmpliconSeq_script.py
       e4b4b52 R-tools - correcting deseq2.R with missing curly bracket
       3ed5904 R-tools : updated help of deseq2.R
       0209640 Python tools - minor formating in AmpliconSeq_script.py
       08591da Python tools - updated parseTrinotateOutput to Python3 syntax
       53e4353 Python tools - correct rrnaBAMconter python3 syntax

2.10.10        Wed Sep 28 12:51:35 2022 -0400        9 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       c316557 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       fd73915 Python tools : CpG_coveragStats now compatible with Python3

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      5 commits

       9204711 Merged in DiffBInd_hotfix (pull request #28)
       7701f42 Debug latest commit
       5595842 Adding contrastnb param for cit test for chipseq - atacseq pipeline
       2303dc7 Reading chrList from input as a string being like "chr1,chr2" (ini file for GenPipes)
       1a84904 Changing awk for allowing genomes to have "-" within its name

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       200a87b added anno_fuse for rnaseq pipeline

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       b49abc9 Merged in rna_seq_update (pull request #27)

2.10.9        Fri Jun 24 18:40:14 2022 +0000        5 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       66770a5 Version bump to 2.10.9
       8a82bca Python tools - Run Processing metrics : improved alignment metrics function

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       7685c96 methylKit.R edited online with Bitbucket

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       81ba9ec added a cutoff to remove based on readcounts

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       aa450a1 Merged in methylkit_diff_methy (pull request #26)

2.10.8        Thu Jun 16 17:00:15 2022 +0000        3 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       7e2580b Version bump to 2.10.8
       c3b5959 Python Tools : run processing metrics improved
       c969134 Version bump to 2.10.7

2.10.7        Tue Jun 7 17:26:05 2022 -0400        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       7cc985c Run Procesisng : setting verifyBamID output as non-mandatory

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      1 commits

       af6e5b1 Version bump to 2.10.6

2.10.6        Fri Jun 3 19:17:57 2022 +0000        3 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       7dff0c5 Version bump to 2.10.5

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       ff40611 Merged in DiffBInd_hotfix (pull request #25)
       8705a96 Debugging asva.R

2.10.5        Fri May 6 12:21:07 2022 -0400        8 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      6 commits

       8b44038 Python Tools - RunRprocessingMetrics2json : allowing absence of verifyBamID metrics
       29a69bf Version bump to 2.10.5
       90cc30e Python Tools - updated preprocess.py to avoid empty lines insertion in the header
       716a3a1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       f63a433 Python Tools - updated preprocess.py for Python3 compatibility
       81b8095 Version bump to 2.10.4

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       a11f9cf Merged in verifybamidfix (pull request #24)

  Rob Syme <rob.syme@gmail.com>      1 commits

       91089e2 Accept both tsv and selfSM extension for VerifyBamID

2.10.4        Mon May 2 16:22:17 2022 +0000        3 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       13eb97d Version bump to 2.10.3

  Édouard Henrion <henrione@narval2.narval.calcul.quebec>      2 commits

       9a0c243 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       890ef0b Python tools - updated fixVS2VCF.py and format2pcgr.py to be python3 compatible

2.10.3        Mon Apr 11 14:48:08 2022 -0400        4 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       6387aa9 Version bump to 2.10.2

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      1 commits

       eae38f3 Merged in DiffBInd_hotfix (pull request #23)

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       cdfd52a Added FDR and P-value cutoff to all the functions

  Pubudu Nawarathna Mudiyanselage <pubudu@narval1.narval.calcul.quebec>      1 commits

       99dfd96 Removed comments

2.10.2        Tue Apr 5 12:37:33 2022 -0400        6 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       9a13ef3 Version bump to 2.10.1

  Édouard Henrion <henrione@narval2.narval.calcul.quebec>      1 commits

       23be1ae PYTHON TOOLS - Updating fixVS2VCF.py to run with Python3

  Édouard Henrion <henrione@narval3.narval.calcul.quebec>      1 commits

       169f9d4 PYTHON TOOLS : updating dict2BEDs.py so it's compatible to both Python2 and Python3

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       02f9aec rrnaBAMcounter.py edited online with Bitbucket

  Robert Syme <rob.syme@gmail.com>      1 commits

       ffc257e Merged in rrnaBAMCounter_pythonupdate (pull request #22)

  Rob Syme <rob.syme@gmail.com>      1 commits

       6569fa5 Remove unused imports and convert print to python3

2.10.1        Wed Mar 9 15:49:27 2022 -0500        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       2f2b163 Python Tools - Run Processing : updating runProcessingMetrics2Json.py with same format as pipeline
       8b349f5 Version bump to 2.10.0

2.10.0        Thu Mar 3 15:15:38 2022 -0500        7 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       d8fff50 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       6af9e78 PYTHON TOOLS - Run Procssing : updates and fixes in runProcessingMetricsToJson.py + removing useless run processing tools
       540404c Python tools - Run Processing : adding tool to push metrics of run processing into a standardized report in SON format
       d9145c7 Version bump to 2.9.1

  Paul Stretenowich <paul.stretenowich@gmail.com>      1 commits

       c731ca7 transfer.sh - Testing removing the mkdir before rsync

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       9beef29 Merged in data_transfer_script (pull request #21)
       3ca743d Updating transfer tool for syncing symlinks

2.9.1        Fri Dec 17 10:47:28 2021 -0500        3 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       bf75bd7 Version bump to 2.9.0

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       8687a39 fixed a typo
       08641a3 added pca plots to the report improved to analyse ATAC-seq peaks

2.9.0        Tue Dec 7 11:06:54 2021 -0500        2 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       b4d1fa5 MUGQIC_TOOLS - adding script to combine metrics from fgbio DemuxFastqs
       605ae14 Version bump to 2.8.2

2.8.2        Tue Nov 30 11:22:49 2021 -0500        6 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       6b7b175 Updated VERSION prior new release
       1f88a55 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       a8ce1bd PYTHON TOOLS - Run Processing : improved metrics in runProcessingSampleReport.py
       43b17d9 Version bump to 2.8.2-beta
       e01530f Version bump to 2.8.1

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       b5cb280 adding signal_to_noise in epiQC pipeline

2.8.1        Wed Oct 27 16:17:08 2021 -0400        6 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       e61de6e Version bump to 2.8.1-beta
       5e1c9d6 Version bump to 2.8.0

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       7f606de MGUQIC_TOOLS - ampliconSeq Tools : updated use of tmp in cleanOTUtable.sh

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       39f2955 modified chipseq DiffBind report Fixed text decorations
       eb08eaf minor fix to the path name
       cba6fda Added KOPathawayAnalysis.R in the seq2fun protocol of the de-novo RNA-seq pipeline

2.8.0        Wed Oct 13 13:17:49 2021 -0400        11 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       8752483 Version bump to 2.7.1-beta
       cf1fd45 Version bump to 2.7.0

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      3 commits

       3df3cf9 Fix to the rsync copy issue provoked by the slash at the end of the input arguemnt (fix by Paul S.)
       ae6effa Fixed but so that script errors out correctly when no arguments are provided
       bb4e920 Fixed permissions issue

  jgalvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       3d972d0 Added first version of data transfer script tools/transfer.sh

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       1d65a64 Merged in data_transfer_script (pull request #20)
       82f6eed Merged master into data_transfer_script

  P-O Quirion <pioliqui@gmail.com>      1 commits

       cc84f72 removing extra \

  Robert Syme <rob.syme@gmail.com>      1 commits

       a4916fe Merged in goseq_fix (pull request #19)

  Rob Syme <rob.syme@gmail.com>      1 commits

       7a21f95 Fixes a bug caused by reading DGE results with NAs

2.7.0        Wed Aug 18 14:14:48 2021 -0400        9 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       3ff4d3d updating version for release

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      2 commits

       fdf9d17 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       3fa0c51 Version bump to 2.6.1

  Paul Stretenowich <paul.stretenowich@gmail.com>      2 commits

       85ff693 Fixing file name finding
       82aae03 Changing name

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       f71165c Merged in IHEC_chipseq_metrics_fix (pull request #18)
       12a2e94 Merged master into IHEC_chipseq_metrics_fix

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       4b13b08 cleanOTUtable.sh edited online with Bitbucket

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       f6829c2 added epiqc.py file for epiQC report generation and signal_noise.py in the epiqc pipeline

2.6.1        Mon Jul 26 19:22:14 2021 +0000        2 commits

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      1 commits

       63e3b93 sequenza fixes

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       6c74997 Merged in tumor_pair (pull request #17)

2.6.0        Fri Jun 18 17:44:38 2021 +0000        8 commits

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      1 commits

       0a770fa Version bump to 2.6.0

  ehenrion <edouard.henrion@mcgill.ca>      2 commits

       5ea76f9 MUGQIC_TOOLS -  bash tools : fixes Issue #1
       dd2eee5 MUGQIC_TOOLS - bash tool : added 'ls' command in cleanOTUtable.sh to ensure 'copy' is over prior 'moving'

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       b5580ca Removed the second if ladder
       92da41b added paramter to change the comparison method. e.g deseq2 or edger
       cc5a254 Decorated html output, center aligned plots, added figure captions and named chunks
       dba8c2b Added some comments explaining the code
       b407204 Added DiffBind.R file used in chipseq differential binding

2.5.0        Tue Mar 30 16:24:40 2021 -0400        155 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       e30895d Version bump to 2.4.0

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      3 commits

       a613bdd Version bump to 2.5.0
       171104d Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       424214c Version bump to 2.5.0

  Paul Stretenowich <paul.stretenowich@gmail.com>      148 commits

       8550435 Fixing typo
       db28643 Fixing typo
       f9edca2 Fixing typo
       d7c0d1c Cleaning and fixing tss script
       3450864 Merge branch 'IHEC_chipseq_metrics_fix' of https://bitbucket.org/mugqic/mugqic_tools into IHEC_chipseq_metrics_fix
       d27f908 Fixing IHEC chipseq metrics
       affddba Fixing IHEC chipseq metrics
       92054ae Fixing IHEC chipseq metrics
       68ebf59 Fixing IHEC chipseq metrics
       8984a6a Fixing IHEC chipseq metrics
       bc7b86f Fixing IHEC chipseq metrics
       1c08b4e Fixing IHEC chipseq metrics
       43fdaf5 Fixing IHEC chipseq metrics
       5ef2251 Switching from design to readset
       f7bf2a9 Switching from design to readset
       d3f5cf8 Switching from design to readset
       75e42c2 Fixing IHEC chipseq metrics
       accd9d5 Fixing IHEC chipseq metrics
       931ee37 Fixing IHEC chipseq metrics
       b37f0a7 Fixing IHEC chipseq metrics
       db2335c Fixing IHEC chipseq metrics
       9d5dcc7 Fixing IHEC chipseq metrics
       7486d44 Fixing IHEC chipseq metrics
       13b932a Merge branch 'IHEC_chipseq_metrics_fix' of https://bitbucket.org/mugqic/mugqic_tools into IHEC_chipseq_metrics_fix
       27126a0 Fixing IHEC chipseq metrics
       71e62a2 Fixing IHEC chipseq metrics
       72d7c2d Fixing IHEC chipseq metrics
       e4b1e4f Fixing IHEC chipseq metrics
       22ce29f Fixing IHEC chipseq metrics
       1c536e4 Fixing IHEC chipseq metrics
       5b1df8d Fixing IHEC chipseq metrics
       effffc2 Fixing IHEC chipseq metrics
       10772dd Fixing IHEC chipseq metrics
       35287e3 Fixing IHEC chipseq metrics
       8b9788c Fixing IHEC chipseq metrics
       a84d4b9 Fixing IHEC chipseq metrics
       bb103de Fixing IHEC chipseq metrics
       231e2cc Fixing IHEC chipseq metrics
       193906a Fixing IHEC chipseq metrics
       2cd3a95 Fixing IHEC chipseq metrics
       7830440 Fixing IHEC chipseq metrics
       3c40592 Fixing IHEC chipseq metrics
       5fa6b2f Fixing IHEC chipseq metrics
       4ed3980 Fixing IHEC chipseq metrics
       223592e Debug
       702a5f9 Switching from design to readset
       62a4bba Switching from design to readset
       da852e1 Switching from design to readset
       be76982 Fixing IHEC chipseq metrics
       8908e84 Switching from design to readser
       e6678d2 Switching from design to readser
       1afb177 Switching from design to readser
       c4aa908 Switching from design to readser
       231ba82 Switching from design to readser
       475eeb1 Switching from design to readser
       693109d Switching from design to readser
       69abdcd Switching from design to readser
       287ed42 Switching from design to readser
       7faa335 Switching from design to readser
       fee2ab9 Switching from design to readser
       adbdd4b Switching from design to readser
       8f0648e Switching from design to readser
       db525ec Switching from design to readser
       7f445c5 Switching from design to readser
       581f387 Switching from design to readser
       c6199ea Switching from design to readser
       dfffe1d Switching from design to readser
       311c4ff Switching from design to readser
       ff14a6c Switching from design to readser
       9b97e77 Switching from design to readser
       a5d777e Switching from design to readser
       d9fb8a7 Switching from design to readser
       b5d2f9d Debug design file change chipseq pipeline
       e0cbfda Fixing IHEC chipseq metrics
       c9c6ec9 Fixing IHEC chipseq metrics
       e69f170 Fixing IHEC chipseq metrics
       f9ca01c Fixing IHEC chipseq metrics
       64235fe Fixing IHEC chipseq metrics
       606d466 Fixing IHEC chipseq metrics
       410c501 Fixing IHEC chipseq metrics
       8afc660 Fixing IHEC chipseq metrics
       249a6b5 Fixing IHEC chipseq metrics
       5ac8c8d Fixing IHEC chipseq metrics
       4127e4e Fixing IHEC chipseq metrics
       bc6879e Debug
       2bd980e Fixing IHEC chipseq metrics
       fd9fbb4 Fixing IHEC chipseq metrics
       36ce67a Fixing IHEC chipseq metrics
       7bc76b0 Fixing IHEC chipseq metrics
       082131c Fixing IHEC chipseq metrics
       053bc5d Fixing IHEC chipseq metrics
       251bd82 Fixing IHEC chipseq metrics
       a71fb4b Fixing IHEC chipseq metrics
       29780fa Fixing IHEC chipseq metrics
       50cd9aa Fixing IHEC chipseq metrics
       42477c9 Fixing IHEC chipseq metrics
       c0e7159 Fixing IHEC chipseq metrics
       e23e896 Fixing IHEC chipseq metrics
       f35b761 Fixing IHEC chipseq metrics
       5929748 Fixing IHEC chipseq metrics
       78da539 Fixing IHEC chipseq metrics
       260750b Fixing IHEC chipseq metrics
       b98693c Fixing IHEC chipseq metrics
       a3087d2 Fixing IHEC chipseq metrics
       086de46 Fixing IHEC chipseq metrics
       affff8a Fixing IHEC chipseq metrics
       b72e3ae Debug
       ac1c9fd Switching from design to readset
       802242c Switching from design to readset
       b5f22a8 Switching from design to readset
       2847309 Fixing IHEC chipseq metrics
       13aacad Switching from design to readser
       15c8ff9 Switching from design to readser
       6617b09 Switching from design to readser
       e22d822 Switching from design to readser
       9d3cd7c Switching from design to readser
       cda9fff Switching from design to readser
       a738c8e Switching from design to readser
       041abc3 Switching from design to readser
       bbd6a4d Switching from design to readser
       e0e96d4 Switching from design to readser
       696a402 Switching from design to readser
       19a380f Switching from design to readser
       6226787 Switching from design to readser
       efdf606 Switching from design to readser
       aa35b4c Switching from design to readser
       712e471 Switching from design to readser
       17bdc05 Switching from design to readser
       29f6f3c Switching from design to readser
       03fcd21 Switching from design to readser
       9206d7a Switching from design to readser
       581600a Switching from design to readser
       02a308b Switching from design to readser
       ab9b34e Switching from design to readser
       bf7a6e5 Debug design file change chipseq pipeline
       e060522 Fixing IHEC chipseq metrics
       ccbb492 Fixing IHEC chipseq metrics
       e1bfb99 Fixing IHEC chipseq metrics
       c6be1a6 Fixing IHEC chipseq metrics
       5442b15 Fixing IHEC chipseq metrics
       9f40a06 Fixing IHEC chipseq metrics
       cd5912f Fixing IHEC chipseq metrics
       b6abcf1 Fixing IHEC chipseq metrics
       644d458 Fixing IHEC chipseq metrics
       e13f7c0 Fixing IHEC chipseq metrics
       66217c5 Fixing IHEC chipseq metrics
       d98d7e0 Debug
       96dd552 Fixing IHEC chipseq metrics

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      3 commits

       22ea87d Merged in IHEC_chipseq_metrics_fix (pull request #15)
       c99130a Merged master into IHEC_chipseq_metrics_fix
       6aaaaa1 Merged master into IHEC_chipseq_metrics_fix

2.4.0        Mon Mar 15 15:27:18 2021 -0400        13 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      8 commits

       6a10cc4 MUGQIC_TOOLS - python tools : updated convertFastqHeaders to convert already demultiplexed fastq files only
       7fa2afd MUGQIC-TOOLS - Python Tools : update convertFastqHeader with use of Bio.SeqIO.QualityIO.FastqGeneralIterator instead of Bio.SeqIO.parse to try speeding things...
       647c5e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       7614df1 pytohn tools - convertFastqHeader - corrected output record handling
       2cc1b36 adding tools for run processing pipeline + fastq header conversion tool
       f6bdd97 adding tools for run processing pipeline + fastq header conversion tool
       304278f Version bump to 2.3.3-beta
       5a3e071 Version bump to 2.3.2

  ehenrion <edouard.henrion@mcgill.ca>      2 commits

       cc43663 Python tools - convertFastqHeaders : corrected typo
       bf2c709 convertFastqHeaders.py edited online with Bitbucket

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       e7f4fdc add fastqIcounter.py from my perso repo

  Robert Syme <rob.syme@gmail.com>      2 commits

       1a983d8 Merged in Robert-Syme/ensure-that-mgi-library-name-parsing-cap-1615835329509 (pull request #16)
       7be259e Ensure that MGI library name parsing captures letter prefix.

2.3.2        Fri Apr 17 13:26:36 2020 -0400        8 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       a6f9432 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       869d24a Version bump to 2.3.1

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      4 commits

       40f63d8 Added a script to run BLAST for Nanopore QC
       6f4c3fb Merge branch 'nanopore_jhg' of bitbucket.org:mugqic/mugqic_tools into nanopore_jhg
       e977d55 Added a tool to trim nanopore reads to a set length of 1000bp
       3ca36d4 Added a tool to trim nanopore reads to a set length of 1000bp

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      2 commits

       7c066b9 Merged in nanopore_jhg (pull request #14)
       368915d Merged master into nanopore_jhg

2.3.1        Mon Mar 30 09:34:18 2020 -0400        4 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       d355b75 Version bump to 2.3.1-beta
       b4e91eb Version bump to 2.3.0

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      2 commits

       945d76e corrected typo in printing new restructured matrix size
       31bdf7a Modified script to take rawRN interaction matrices. Added additional codes to re-structure the code

2.3.0        Fri Feb 21 15:19:26 2020 -0500        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       aab730b correcting file permissions to some python tools
       d9ba47c Version bump to 2.2.5-beta
       5cba859 Version bump to 2.2.4

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      1 commits

       4a2b153 Added a new python script to extract a certain number of bases from long nanopore reads, for use with BLAST QC.

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       8b61def Merged in nanopore_jhg (pull request #12)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      10 commits

       4febc34 Merged in IHEC_chipseq_metrics_fix (pull request #11)
       2b4d612 IHEC_chipseq_metrics_max.sh edited online with Bitbucket - Fix typo
       208edcd IHEC_chipseq_metrics_max.sh edited online with Bitbucket - Ed comment
       b57833b IHEC_chipseq_metrics_max.sh edited online with Bitbucket - Fixing sambamba issue
       e3ad940 Using sambamba to speed up
       6108b26 Using sambamba to speed up
       cfd4925 Fixing metrics
       89d50d5 Fixing typo
       11e5f3d Merged in Paul-Stretenowich/getfastabinedgcpy-edited-online-with-bit-1558131034321 (pull request #10)
       af4f10f Correcting minor bug issue

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      1 commits

       b85eb34 Added hicrep.R to R-tools which does the hicrep analysis

  rami.coles@mail.mcgill.ca <rcoles@abacus3.ferrier.genome.mcgill.ca>      1 commits

       a83ea10 Added tool to compute signal to noise on ChIP-Seq signal tracks and EpiQC report tool

2.2.4        Thu May 16 11:01:29 2019 -0400        13 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       4e1b0f4 added tool to create the cluster preset for SMRTLink
       93ea18e Version bump to 2.2.4-beta
       af0d356 Version bump to 2.2.3

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       9064f33 Merged in fix_asva_plot_error (pull request #9)
       49fa3dd Merged in fix_createBaitMapFile (pull request #8)

  P-O Quirion <pioliqui@gmail.com>      2 commits

       3a6188e fix for the plotting at the end
       2159261 rewrite escaping

  Robert Syme <rob.syme@gmail.com>      5 commits

       e9d2b20 Merged in ballgown_fix (pull request #6)
       8cde138 Explicit use of dplyr::contains
       6380b43 Non-essential stylistic changes
       fc52621 Ensures that design files that use empty cells instead of '0' are read correctly.
       f91478c Ignore basic R project files (just in case).

  Robert Syme <rsyme@seychelles.genome.mcgill.ca>      1 commits

       d47ccdb Remove hard-coded output name and use file.path

2.2.3        Wed May 1 15:46:36 2019 -0400        6 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       96eef34 commit before release
       630bca1 Version bump to 2.2.3-beta
       d13833a Version bump to 2.2.2

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       5859f60 Merged in fix_createRmapFile (pull request #7)

  P-O Quirion <pioliqui@gmail.com>      2 commits

       7654b6f update createRmapFile.sh
       b04d6ef update createRmapFile.sh

2.2.2        Mon Mar 11 15:24:56 2019 -0400        7 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       52bd87e commented some part of code used for testing purposes : were only considering part of input to speed-up the process
       27bee98 MethylSeq DMR analysis - flip the case and control assignment in methylKit.R
       b39724f Version bump to 2.2.2-beta
       23676fa Version bump to 2.2.1

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      2 commits

       76b4d0c Corrected error caused by sleuth version update, as well as other minor mistakes. Script now produces PCA plots as well as heatmap of top 40 genes by FC.
       09b8a1b Corrected tibble error with the sleuth script

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       5fe716e Merged in rnaseq_light_jhg (pull request #5)

2.2.1        Thu Feb 28 11:43:00 2019 -0500        14 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      2 commits

       9ed419b Version bump to 2.2.1-beta
       b06e519 Version bump to 2.2.0

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      7 commits

       83bfc13 Fixed merge issues with sleuth script
       e63ff6a Merge branch 'master' into rnaseq_light_jhg
       04c4d4f Fixed merge conflit with kallisto script
       a6f6d4b Added corrections to the kallisto script, to incorporate logs (for MultiQC) and reordering parameters to avoid conflicts.
       f2bf76e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_tools
       f9e276a Incorporated updates to the kallisto script to save log (for MultiQC) and to fix order of parameters
       2fa54d3 Added initial sleuth scritp (incorporating corrections by Eloi M.)

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      2 commits

       5e740d1 Merged rnaseq_light_jhg into master
       98c7f70 Merged master into sleuth

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       38adb2c Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       32c293b remove bug for OntargetReadsDedup calculation

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       d0b724f add python script to add genotypes in strelka2 vcf

2.2.0        Tue Dec 18 10:42:31 2018 -0500        12 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       c37d23e added the loading of 'caTools' library in run_spp.R - avoid issues when loading spp library
       ceda616 Added methylKit.R
       bb4a410 Version bump to 2.1.13-beta
       3bc4e64 Version bump to 2.1.12

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      4 commits

       828509f Corrected added bug on line 92
       83eb168 Clarified R versions commment
       ac76e73 Added R version to header of script
       dafd91f First commit of ballgown R script, compatible with StringTie protocol of GenPipes RNA-seq pipeline

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      2 commits

       0206e2f Merged ballgown into master
       ccaf25f Merged in ballgown (pull request #2)

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       5bf8dd8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       c03602d correct bugs and add on target stats

2.1.12        Fri Nov 16 15:16:04 2018 -0500        18 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      12 commits

       95d7d73 Version bump to 2.1.13-beta
       0f58aa5 updated CHANGELOG
       1912bd5 updated CHANGELOG
       aa27460 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       20b7ee0 Adding R-tools for Sleuth & FusionMetaCaller
       03ef687 Version bump to 2.1.13-beta
       a23b5da Version bump to 2.1.12
       5086606 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       581ebf7 Adding asva.R - used by of ampliconseq dada2 protocol
       2999f1f added ampliconLengthFile parameter to asva.R
       0d53986 Ampliocon-Seq pipeline - dada2 protocol : add asva.R script to R-tools
       a7c1ef0 Version bump to 2.1.11

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      4 commits

       13f8106 corrected typo in R-tools/asva.R (missing parenthesis) and updated perl-tools/methylProfile.bismark.pl
       83a5fb5 Corrected missing parenthesis in asva.R
       25cce0d add asva.R : used by the dada2 protocol of the ampliconseq pipeline
       0d11c07 minor corrections for the MethylSeq coverage stats

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      1 commits

       15c2429 added clean_otu.sh script to mugqic_tools, used in ampliconseq pipeline to clean the OTU tables (remove lines containing undesired characters e.g. division, OP3, WS6...)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       bbfa8a2 update IHEC_metrics for methylSeq for targeted

2.1.11        Mon Mar 26 14:22:46 2018 -0400        20 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      9 commits

       8bcae89 changed permission to 'a+x' for some R-tools used by smallRNA pipeline
       9c54695 added new bash tools used by HiC-Seq pipeline
       fab5e28 minor updates regarding shebang of our bash tools and updated indentation
       44fa837 minor indentations updats + removed flagstat calls from IHEC_methylseq_metrics since they are now made earlier in the methylseq pipeline
       e2ad17b Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       7e098c3 smallRNA pipeline - added some more R tools to handle data for/from mirdeep2
       56f6f34 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       24d2b3a smallRNA pipeline - added R tool to prepare the inputs for mirdeep2
       4f7b3d0 updates brought to version 2.1.10 dump

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      8 commits

       72ab8a2 in mergeKnownAndNovelMiRNA.R: change variables name cleanup unused commands
       e0cf384 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       dd065bd in plotSpikeInCount.R: add a column for percentage among tags and change column names of summary file
       7310bf0 in abundanceTranscript2geneLevel.R: add parameter to ignore taxon version which fix error when gtf and fasta use different annotations
       691ac26 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       2acd2b6 in estimateSpikeInCount.sh, fix input file name when running plotSpikeInCount.R
       802cdf9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b2b91a7 in estimateSpikeInCount.sh, filterBlastOutput.py and plotSpikeInCount.R: add .tsv file extension

  Rola Dali <rola.dali@mail.mcgill.ca>      3 commits

       74734b5 updating RobusTAD for TAD calling
       a54b860 adding RobusTAD for hicseq.py
       8c907d1 adding MT reads to IHEC metrics

2.1.10        Thu Dec 21 15:49:09 2017 -0500        8 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       970f924 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       16c8203 dict2bed - correction to avoid error message : "global name 'args' is not defined"
       791827b Version bump to 2.1.10
       89605c5 BFXDEV-674 - MethylSeq pipeline - updated metrics scprits for MethylSeq pipeline : standardized header and some reviewed calculations
       f421c9b Version 2.1.9 updated

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      2 commits

       e433367 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       a1719b9 in plotSpikeCount.R remove unused variable pdf_out

  mmichaud <marc.michaud@mail.mcgill.ca>      1 commits

       4cf22bd Add CountIlluminaBarcodes in java-tools; used in illumina run processing pipeline.

2.1.9        Thu Dec 14 11:27:05 2017 -0500        91 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      34 commits

       ac44095 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       68e8be3 very minor updates on the header of the methylseq metrics report
       ad0c35d BFXDEV-674 - updated IHEC metrics report script : generalized way to find bismark alignment reports when computing the human conversion rate
       3dac5d2 Version bump to 2.1.9
       fc8753d minor updates in metrics (& IHEC metrics) scripts for MethylSeq pipeline
       7558a25 new CHANGELOG.md
       5b25386 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       496f8fb BFXDEV-674 - review the metrics step, especially the IHEC metrics step, and updated header labels to be more standard
       21b5428 new CHANGELOG.md
       a8c8180 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       e9bb399 updates on IHEC metrics for ChipSeq
       fd6b95b new CHANGELOG.md
       f68d795 added run_spp.R to make proper use of 'spp' R library
       88e1287 new CHANGELOG.md
       85f5cd4 adding tools/IHEC_chipseq_metrics_max.sh
       1ff0790 new CHANGELOG.md
       c1ed44a Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       09fc149 BFXDEV-674 - updated MethylSeq metrics scripts to reflect changes in MethylSeq pipeline regarding the deduplication step (which now uses Picard MarkDuplicates)
       a79b874 new CHANGELOG.md
       f27476c BFXDEV-674 - added a comment mark to the header line of the output combined file
       2c1edba BFXDEV-674 - corrected unset variable in cpgStats.sh
       b49cd5d BFXDEV-674 - added tools/cpgStats.sh for CpG metrics computing within the methylSeq pipeline
       ce53fe2 BFXDEV-674 - adding bash tools for methylSeq metrics
       019df95 releasing README and CHANGELOG for verion 2.1.9
       a119ec4 added header line in GC content file
       957de36 BFXDEV-674 - perl tool methylProfile.bismark.pl minor update
       bf43c15 while working on BFXDEV-172, did a minor review or goseq.R and edger.R : updated indentation for more clarity when reading the code
       1da92ed BFXDEV-674 - added python tool to assess the CpG coverage stats in the MethylSeq pipeline
       f40eaf3 BFXDEV-674 - added R tool for the calculation of the GC bias in the MethylSeq pipeline
       2e271d4 BFXDEV-674 - added python tool to calculate GC content from a given fasta file (e.g. reference), mainly used to produced GC content reference file needed by the MethylSeq pipeline when calculating the GC bias
       e2f3ed7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b7a09d5 MethylSeq - added perl tool to compute methylation profile
       e8689b8 BFXDEV-172 - implemented DESeq2 R-tool and added it to mugqic_tools
       b1efd80 updated CHANGELOG.md from 2.1.8 release

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      20 commits

       0afad0a in plotSpikeCount.R, round identity percentage for clarity
       b8f1cfa in plotSpikeCount.R, remove border of barplots
       f489b1c in plotSpikeCount.R, move percentaage axis to the right, add grid to plots
       035f999 In plotSpikeInCount, add identity barplot and grou tag count and tag percentage into one plot
       ce65309 In plotSpikeInCount, add option beside=T to barplots
       cb0f8cd add plot showing distribution of spike in tag
       5a57302 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       3a5a7f7 fix an error with query start and end when filtering blast output
       06beb03 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       2938dcf adding scripts for spike-in step in illumina pipeline: plotSpikeInCount.R, filterBastOutput.py, estimateSpikeInCount.sh
       5892d79 remove a call to a non defined variable
       9e6f5c4 fix unused parameter in abundanceTrancript2gene.R
       6040b75 added new parameters to kallisto script to handle single reads
       54a0f98 replaced mentions of samples to readsets
       717ffcc fix spelling
       5638e15 fix output of abundanceTranscript2geneLevel function
       2a7446d added direct call to R_TOOLS variable
       cbf698f added a step for merging individual abundance files
       1d393c1 Add R function to convert transcript to gene level abundance
       aa7ae30 added new function for RNAseq_light

  eloi.mercier@mcgill.ca <eloi.mercier@mcgill.ca>      1 commits

       9de32e1 Merged in RNAseq_light_dev (pull request #1)

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       b84d9e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       b8acd7b IHEC_chips metrics tool - correct typo
       8591e6f tool -  ihec metrics debug -  BFXDEV-675
       a25d589 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       2a359af  ChIPseq-   debug  IHEC metrics tools - BFXDEV-675

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       d76a8f8 tools - cpgStats add more info to stdout and correct for the exit code 1 when grep do not find a match
       ca5603b tools - Methylseq ihec metrics adjustement to the new dedup step
       a21eb00 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       5010cb8 tools - Methylseq ihec metrics adjustement to the new dedup step
       e5d11d8 add tool to recoded sample name in the reedset file based on project.nanuq file (library Name, lane, run)
       2631f72 change ihec methylseq metrics to work on single sample and add GC_bias
       ffa02e5 python tools - add script to extract cancer SNVs from gemini db to xls
       d8cc34b  RNAseq & ChIPseq-   Update (chip) and debug (Rna) IHEC metrics tools - BFXDEV-668  - BFXDEV-675
       456de6c tools -  implement IHEC ChIPseq metrics script - BFXDEV-675
       3fb02cc tools -  implement IHEC RNA metrics script - BFXDEV-668
       3d2e2af Add IHEC rnaseq metrics generation script
       3a8ab8d update axiom dev suite
       a1ea573 correct plateQC graphs resolution
       787a517 correct plateQC graphs

  Mathieu Bourgey <mbourgey@cedar5.cedar.computecanada.ca>      1 commits

       4493877 cadjust permission of tools/IHEC_chipseq_metrics.sh tools/IHEC_rnaseq_metrics.sh - BFXDEV-668  - BFXDEV-675

  mmichaud <marc.michaud@mail.mcgill.ca>      1 commits

       9df614a Allow a minimum subsampling threshold (1x10e-6) to blast even with a lot of source reads. BFXDEV-698

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      1 commits

       aaf5581 Update sequenza script and add sCNAphase script

  Rola Dali <rola.dali@mail.mcgill.ca>      14 commits

       6069ee9 changing ihec_metrics to use chip_type N/B instead of chip mark
       a868c97 IHEC_chipseq_metrics_max.sh edited online with Bitbucket
       7bb9b0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       42ef372 fixing rRNA matrics for IHEC
       5ccb9d2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       0199656 IHEC format changes
       a98798a chipseq modification: add rawreads
       73fbca8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       e66dc77 ihec_metrics_rnaseq.py for rnaseq metrics
       987da9e rosolving merge conflitct
       28d96a9 adding NSC and RSC
       092451b Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       c1040bf adding CreateHicFileInput.sh for .hic file creation in hicseq BFXDEV-670
       aa6c18f added hicseq related scripts BFXDEV-670

2.1.8        Mon Apr 24 09:55:28 2017 -0400        18 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       edfb29a commiting the CHANGELOG before release

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      10 commits

       9ec9b5f correct ggplot
       e638845 correct ggplot
       eafb4fc correct typo
       7f025be correct typo
       1b9f6b8 correct typo
       9f98c73 correct typo
       dd2ba39 correct typo
       170f2b4 correct typo
       5e77a35 correct typo
       9ccb828 correct module incompatibility

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       70af3b5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       69cce81 update gene titan
       068bd1b Imporve Genetitan processing script
       4c049c5 Imporve Genetitan processing script
       446f05d resolve conflict merging
       94240fc Updat Axiom tools for more QC
       3f461c4 add somatic signature R tools

2.1.7        Fri Jan 20 16:31:19 2017 -0500        19 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      10 commits

       2689de8 MUGQIC_TOOLS - added fixVS2VCF.py to python-tools as it is needed in the tumor_pair pipeline
       eae6a16 commit prior to mugqic_tools 2.1.7 release
       aa281e5 deseq.R - added '-l => local fit' parameter explanation to the tool manual
       79f0a77 Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       91e39ac DESeq.R - using 'locfit' to allow locat fit instead of just having to use parametric dispersion
       0db2adb Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       1791395 DESeq & Edger - explicitely importing 'methods' library because Rscript does not
       5f65c6e Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       0e76111 metaMArkerSeq - correct shebang for heatmap Rscript creation
       7f27678 ampliconSeq - krona debug

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       6f9a78d Add R script to filter Gene Titan processed data
       a21249f Add bash script to process gene titan
       1e3420b Merge branch 'master' of bitbucket.org:mugqic/mugqic_tools
       ff2f40d R-tools - remove typos in RunSequenza_analysis.R
       18d5720 make runBlast to support no empty fastq 3rd line
       822983a adding sequenza analysis script

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      1 commits

       efc56e0 Adding vcf2bed.pl script for varscan2 of tumor_pair pipeline

  Tushar Dubey <tdiitg@gmail.com>      2 commits

       026180c removed subset line from script PopSV
       8e0170c Added PopSV script(beta)

2.1.6        Tue May 3 14:55:43 2016 -0400        5 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       a929638 new CHANGE LOG prior to new release 2.1.6
       39ac81a resolve conflicts
       362d40a commiting changes prior to new release
       1c82601 new CHANGELOG.md

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      1 commits

       df62d30 dict2bed - correct bug and add feature - BFXDEV-521

2.1.5        Tue Jan 26 11:57:08 2016 -0500        8 commits

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

