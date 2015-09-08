#!/usr/bin/env python

# -*- coding: utf8 -*-

import re, sys
import operator
import collections
import numpy
import argparse

def bfilter(inputF,outputF,num_species,Tqcovs):

	VinputF = open(inputF,'r')	# Input
	VoutputF = open(outputF,'w')	# Output
	
	dic_taxon={}	# Taxon dictionary
	dic_pident={}	# Identity percentage dictionary
	dic_qcovs={}	# Query coverage dictionary
		
	lines = VinputF.readlines()
	i=0
	
	while i<len(lines):
	
		word = re.split("[\r\t]",lines[i])	#1	
		w_taxon = word[7]
		w_pident = word[5]
		w_qcovs = word[8]
		
		if w_taxon != 'N/A':
		
			if float(w_qcovs)>float(Tqcovs):
			
				word2 = re.split(" ",w_taxon)
				w_genus = word2[0]
				w_species = word2[1]			
				
				if w_genus+'_'+w_species in dic_taxon.keys():
					dic_taxon[w_genus+'_'+w_species]+=1
				else:
					dic_taxon[w_genus+'_'+w_species]=1
					
				if w_genus+'_'+w_species in dic_pident.keys():
					dic_pident[w_genus+'_'+w_species].append(float(w_pident))
				else:
					dic_pident[w_genus+'_'+w_species]=[float(w_pident)]
			
				if w_genus+'_'+w_species in dic_qcovs.keys():
					dic_qcovs[w_genus+'_'+w_species].append(float(w_qcovs))
				else:
					dic_qcovs[w_genus+'_'+w_species]=[float(w_qcovs)]

		i+=1			
	
	dic_taxon_format = collections.Counter(dic_taxon)
	
	VoutputF.write('Species'+'\t'+'Hit Count'+'\t'+'Average Percent Identity'+'\t'+'Average Query Cover (%)'+'\n')
	for hit in dic_taxon_format.most_common(int(num_species)):	# Top taxon
		VoutputF.write(hit[0]+'\t'+str(hit[1])+'\t'+str(int(round(numpy.mean(dic_pident[hit[0]]))))+'\t'+str(int(round(numpy.mean(dic_qcovs[hit[0]]))))+'\n')

	VinputF.close()
	VoutputF.close()
	

def main(argv):
	
	mod=[]
	mod.append('\n%(prog)s -i <blast_result> -s <num_species> -q <min_query_coverage> -o <output>')
	
	parser = argparse.ArgumentParser(prog = 'blast_filter.py',
                                 usage = "\n".join(mod))
                                 
	parser.add_argument('-i', action='store', dest='input_value',
	                    help='Input (blast result).')
	                    
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output.')

	parser.add_argument('-s', action='store', dest='num_species_value',
	                    help='Number of species to keep.')
	                    
	parser.add_argument('-q', action='store', dest='qcovs_value',
	                    help='Query coverage minimum value.')
	                    	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()
	
	# Blastn filter
	
	if results.input_value and results.output_value and results.num_species_value and results.qcovs_value:
		bfilter(results.input_value,results.output_value,results.num_species_value,results.qcovs_value)
	else:
		print 'See help for options: blast_filter.py -h'
									
if __name__ == "__main__":
	main(sys.argv[1:])
	
	