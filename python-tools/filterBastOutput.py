#!/usr/bin/env python
'''
created on February 24, 2017
'''
from optparse import OptionParser
import os,sys
import logging
import re

usage = """ %prog --input_file file.blastres --tag_file tag_database --output_dir output_dir [--perc_identity_min perc_identity_min] [--align_length_min align_length_min] [--mismatches_max mismatches_max] [--skip_blast] [--logging_level log_level --logging_file log_file]
This script parses the output file of blast and return the alignments that satisfy the parameters. """

def main():
    ###################
    #initialization
    ###################
    parser = prep_opt_parser()
    (options,args) = parser.parse_args()
    inputs = verify_inputs(options)
    basename=os.path.basename(options.input_file)

    ###################
    #loggings
    ###################
    if options.log_file:
        log_file=options.log_file
    else:
        log_file=options.output_dir + "/" + basename+".log"
    print("Writing logs into %s" %log_file)
    numeric_level = getattr(logging, inputs.log_level.upper(), None)
    logging.basicConfig(filename=log_file,format='%(levelname)s:%(message)s',level=numeric_level, filemode="w") #print to file

    ###################
    #output files
    ###################
    blast_file=options.input_file
    logging.debug("Blast file %s" %blast_file)
    output_filter_file=blast_file + "." + str(options.align_length_min) + "bp_" + str(options.mismatches_max) + "MM_" + re.sub("\.", "",str(options.perc_identity_min)) + "id.tsv"
    logging.debug("Outout filter file %s" %output_filter_file)

    ###################
    #FILTER
    ###################
    run_blastres_filter(options, blast_file, output_filter_file)

    logging.info("Done")
    return 0


#ensure all input requirements met
def verify_inputs(options):
    if not (options.output_dir):
        logging.critical("Missing output file. Please use -h or --help for input requirements")
        sys.exit(-1)
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if not (options.input_file):
        logging.critical("Missing input file. Please use -h or --help for input requirements")
        sys.exit(-1)
    if not os.path.isfile(options.input_file):
        logging.critical(options.input_file + " is not a valid file.")
        sys.exit(-1)
    if not (options.tag_file):
        logging.critical("Missing tag file. Please use -h or --help for input requirements")
        sys.exit(-1)
    if not os.path.isfile(options.tag_file):
        logging.critical(options.tag_file + " is not a valid file.")
        sys.exit(-1)
    if not (0 <= options.perc_identity_min <= 1):
        logging.critical("perc_identity_min must be between 0 and 1. Please use -h or --help for input requirements")
        sys.exit(-1)
    return options

#prepare parser options
def prep_opt_parser():
    parser = OptionParser()
    #inout files
    parser.add_option("--input_file",dest="input_file",default=None,help="location of the fastq or fastq.gz file to be processed [REQUIRED]")
    parser.add_option("--tag_file",dest="tag_file",default=None,help="location of the file containing the tags, must be a fasta file [REQUIRED]")
    parser.add_option("--output_dir",dest="output_dir",default=None,help="location of the output directory [REQUIRED]")
    #arguments
    parser.add_option("--perc_identity_min",dest="perc_identity_min",type="float", default=0,help="the minimum percentage of identity for the aligments to return (default 0) [OPTIONAL]")
    parser.add_option("--align_length_min",dest="align_length_min", type="int", default=1, help="the minimum length of the alignments to return (default 1) [OPTIONAL]")
    parser.add_option("--mismatches_max",dest="mismatches_max", type="int", default=100,help="the maximum number of mismatches allowed (default 100) [OPTIONAL]")
    #options
    parser.add_option("--logging_level",dest="log_level",default="info", help="set logging level (INFO or DEBUG) [OPTIONAL]")
    parser.add_option("--logging_file",dest="log_file",default=None, help="log file (default input_file.log) [OPTIONAL]")
    return parser

#print subprocess output to logs
def log_subprocess_output(process, pipe):
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        logging.info('%s: %r' %(process,line))

def get_clean_tag_name(tag_name):
    tag_name_no_bracket = re.sub('>', '', tag_name)
    tag_name_clean=tag_name_no_bracket.split("|")[-1]
    return tag_name_clean

#try delete file
def silent_delete_file(file_to_delete):
    try:
        logging.info("Deleting %s" %file_to_delete)
        os.remove(file_to_delete)
    except OSError:
        pass

#run filtering on blast file
def run_blastres_filter(options, blast_file, output_filter_file):
    if not os.path.isfile(blast_file):
        logging.critical(blast_file + " doesn't exist. Verify that you ran BLAST.")
        print("CRITICAL:"+ blast_file + " doesn't exist. Verify that you ran BLAST.")
        sys.exit(-1)

    logging.info("Starting filtering blast file with:")
    logging.info("\tInput file: %s" % blast_file)
    logging.info("\tMinimum alignment length: %i" % options.align_length_min)
    logging.info("\tMinimum identity: %f" % options.perc_identity_min)
    logging.info("\tMismatch max: %i" % options.mismatches_max)

    output_filter_feed=open(output_filter_file,"w")
    row_cpt=0
    tag_found_row=-999
    aligment_pass=False
    all_info_collected=[False,False] #in this order: tag name //, info, query, subject
    count_tags=0

    #collect tag length
    tag_names=[]
    tag_length=[]
    tag_data=open(options.tag_file)
    for tag_line in tag_data:
        tag_entry=tag_line.strip()
        if tag_entry.startswith(">"):
            tag_names.append(get_clean_tag_name(tag_entry))
        else:
            tag_length.append(len(tag_entry))
    tag_dic=dict(zip(tag_names, tag_length))
    tag_data.close()

    blast_alignments=open(blast_file)
    header="Tag\tTag_length\tRead\tRead_length\tAlignment\tIdentity\tMismatches\tQuery_start\tQuery_end\n"
    output_filter_feed.write(header)

    read_name=""
    count_tags=0
    for line in blast_alignments:
        entry=line.strip().split() #qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore

        if (entry[0] != read_name): #compare current name to name of last alignment
            new_read=True
        read_name=entry[0]
        read_length=entry[1]
        tag_name=get_clean_tag_name(entry[2])
        tag_length=entry[3]
        identity=float(entry[4])
        align_size=int(entry[5])
        mismatches=int(entry[6])
        query_start=int(entry[7])
        query_end=int(entry[8])

        if(int(align_size)>=int(options.align_length_min) and float(identity)>=float(options.perc_identity_min) and int(mismatches)<=int(options.mismatches_max)):
            logging.debug("WRITING OUTPUT")
            logging.debug("PASS CRITERIA: %r" % aligment_pass)
            logging.debug("\tAlignment size: %i (required %i)" % (align_size,options.align_length_min))
            logging.debug("\tIdentity: %f (required %f)" % (identity,options.perc_identity_min))
            logging.debug("\tMismatch: %i (maximun %i)" % (mismatches,options.mismatches_max))
            output=tag_name+"\t"+str(tag_length)+"\t"+read_name+"\t"+str(read_length)+"\t"+str(align_size)+"\t"+str(identity)+"\t"+str(mismatches)+"\t"+str(query_start)+"\t"+str(query_end)+"\n"
            output_filter_feed.write(output)
            if (new_read):
                count_tags+=1

        # row_cpt+=1

        # entry=line.strip().split()
        # if (not line.startswith(("Query=",">"," Identities"))): #increase speed by ignoring non relevant rows
        #     continue

        # if line.startswith("Query="): #get the name of the current read
        #     read_name=entry[1]
        #     continue

        # if line.startswith(">"): #grab the tag name
        #     logging.debug(read_name)
        #     logging.debug(line.strip())
        #     tag_found_row=row_cpt
        #     tag_name=get_clean_tag_name(entry[0])
        #     all_info_collected[0]=True
        #     continue

        # if (row_cpt==(tag_found_row+4)):  #4 row later, get aligment info
        #     logging.debug(line.strip())
        #     identity_info=entry[2].strip().split('/')
        #     matches=int(identity_info[0])
        #     align_size=int(identity_info[1])
        #     identity=matches/float(align_size)
        #     mismatches=align_size-matches
        #     if tag_name not in tag_dic:
        #         print("%s no in tag database!" %tag_name)
        #         sys.exit(-1)
        #     missing_bp=tag_dic.get(tag_name)-align_size
        #     # true_mismatches=max(mismatches, mismatches+missing_bp) #58/60bp aligned, 0mm reported => 2 *true* mismatches; also covers tag:58bp, align 58/60bp with 2mm
        #     true_mismatches=mismatches
        #     if(int(align_size)>=int(options.align_length_min) and float(identity)>=float(options.perc_identity_min) and int(true_mismatches)<=int(options.mismatches_max)):
        #         aligment_pass=True
        #     else:
        #         aligment_pass=False
        #     logging.debug("PASS CRITERIA: %r" % aligment_pass)
        #     logging.debug("\tAlignment size: %i (required %i)" % (align_size,options.align_length_min))
        #     logging.debug("\tIdentity: %f (required %f)" % (identity,options.perc_identity_min))
        #     logging.debug("\tMismatch: %i (maximun %i)" % (true_mismatches,options.mismatches_max))
        #     all_info_collected[1]=True

        # if (all_info_collected and aligment_pass):
        #     logging.debug("WRITING OUTPUT")
        #     output=tag_name+"\t"+read_name+"\t"+str(align_size)+"\t"+str(identity)+"\t"+str(true_mismatches)+"\n"
        #     output_filter_feed.write(output)
        #     #reiniziatlize
        #     all_info_collected=[False,False]
        #     tag_found_row=-999
        #     aligment_pass=False
        #     count_tags+=1

    logging.info("Filtering completed. %i alignments matching the parameters found." % count_tags)
    blast_alignments.close()
    output_filter_feed.close()
    return 0

#main
if __name__ == "__main__":
    main_status = main()
    if (main_status==0):
        print "Filtered file successfully generated"
    else:
        logging.critical("Writing of filtered file failed!")
        print("CRITICAL: Writing of filtered file failed!")
        sys.exit(-1)
