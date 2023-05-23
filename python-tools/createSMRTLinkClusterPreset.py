#!/usr/bin/env python

### Edouard Henrion (2019/05/07) - edouard.henrion@computationalgenomics.ca

import os
import sys
import getopt
import xml.etree.ElementTree as ET


def getarg(argument):
    chunk_mode = str(False)
    tmp_dir = ""
    xml_output = ""

    try: 
        opts, args = getopt.getopt(argument[1:],"hct:o:",['help','chunk','tmp=','output='])
    except:
        usage()
        sys.exit(2)

    if len(opts) == 0 :
        usage()
        sys.exit("Error : No argument given")

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        if opt in ("-c", "--chunk"):
            chunk_mode = str(True)
        if opt in ("-t","--tmp"):
            tmp_dir = arg
        if opt in ("-o","--output"):
            xml_output = arg

    if not os.path.exists(tmp_dir) :
        sys.exit("Error - temporary directory not found :\n" + str(tmp_dir))

    print "Chunk Mode is %s" % str(chunk_mode)
    print "Temporary Directory is %s" % str(tmp_dir)
    print "Output will be written in %s" % str(xml_output)

    return chunk_mode, tmp_dir, xml_output
    
def createXmlPreset(chunk_mode, tmp_dir, xml_output):
    root = ET.Element('pipeline-template-preset')
    root.set('id', 'C3G_Preset')
    options = ET.SubElement(root, 'options')

    dist_mode = ET.SubElement(options, 'option')
    dist_mode.set('id', 'pbsmrtpipe.options.distributed_mode')
    dist_val = ET.SubElement(dist_mode, 'value')
    dist_val.text = "False"

    chunk_flag = ET.SubElement(options, 'option')
    chunk_flag.set('id', 'pbsmrtpipe.options.chunk_mode')
    chunk_val = ET.SubElement(chunk_flag, 'value')
    chunk_val.text = str(chunk_mode)

    tempo_dir = ET.SubElement(options, 'option')
    tempo_dir.set('id', 'pbsmrtpipe.options.tmp_dir')
    tempo_val = ET.SubElement(tempo_dir, 'value')
    tempo_val.text = str(tmp_dir)

    tree = ET.ElementTree(root)
    tree.write(xml_output)

def usage():
    print "USAGE : createSMRTLinkClusterPreset.py [option] -o <output_file>"
    print "  -c,               --chunk              :    activate chunk mode"
    print "  -t <folder_path>, --tmp <folder_path>  :    temporary directory"
    print "  -o <file_path>,   --output <file_path> :    output XML file"
    print "  -h,               --help               :    this help \n"


def main():
    print "\n---------------------------------------------------------------------------------"
    print "createSMRTLinkClusterPreset.py creates a preset XML file needed for SMRTLink to run."
    print "The XML file created contains settings for SMRTLink regarding its cluster usage."
    print "This program was written by Edouard Henrion"
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "----------------------------------------------------------------------------------\n"
    chunk_mode, tmp_dir, xml_output = getarg(sys.argv)
    createXmlPreset(chunk_mode, tmp_dir, xml_output)
    
main()
