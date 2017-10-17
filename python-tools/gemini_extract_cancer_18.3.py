import sys
from gemini import GeminiQuery

database = sys.argv[1]
sample_file = sys.argv[2]

dataFile = open(sample_file, 'r')

header = ['chrom',
          'start',
          'end',
          'ref',
          'alt',
          'rs_ids',
          #'filter',
          'cosmic_ids',
          #'qual',
          #'QD',
          'gene',
          'transcript',
          'num_hom_ref',
          'num_het',
          'num_hom_alt',
          #'families',
          'variant_samples',
          'HET_samples',
          'HOM_ALT_samples',
          'pfam_domain',
          'codon_change',
          'aa_change',
          'impact',
          'impact_severity',
          'type',
          'sub_type',
          'call_rate',
          'aaf',
          'in_esp',
          'aaf_esp_ea',
          'aaf_esp_all',
          'in_1kg',
          'aaf_1kg_amr',
          'aaf_1kg_eur',
          'aaf_1kg_all',
          'in_exac',
          'aaf_exac_all',
          'aaf_adj_exac_all',
          'in_omim',
          'clinvar_sig',
          'clinvar_disease_name',
          'clinvar_dbsource',
          'clinvar_dbsource_id',
          'clinvar_origin',
          'clinvar_dsdb',
          'clinvar_dsdbid',
          'clinvar_disease_acc',
          'clinvar_in_locus_spec_db',
          'clinvar_on_diag_assay',
          'gms_illumina',
          'rmsk',
          'in_cpg_island',
          'in_segdup',
          'is_conserved',
          'in_cse',
          'grc',
          'gms_illumina',
          'gerp_bp_score',
          'gerp_element_pval',
          'recomb_rate',
          'cadd_raw',
          'cadd_scaled',
          'fitcons']

samples = list()

for eachLine in dataFile:
    eachLine = eachLine.rstrip("\n")
    samples.append(str(eachLine))
    DP = "DP." + eachLine
    AF = "AF." + eachLine
    GT = "GT." + eachLine
    header.append(str(DP))
    header.append(str(AF))
    #header.append(str(GT))
    header.append(str(eachLine))

print "\t".join(header)

gq = GeminiQuery(database,include_gt_cols=True)
smp2idx = gq.sample_to_idx
query = "SELECT * FROM variants WHERE in_cpg_island = 0 AND in_segdup = 0 AND in_cse = 0 AND rmsk IS NULL" # AND impact_severity <> 'LOW'"

#gt_filter = "(gt_depths).(*).(>=10).(all)"

gq.run(query,show_variant_samples=True)
#gq.run(query,show_variant_samples=True,gt_filter=gt_filter)

for row in gq:
    #print row
    #try:
    #str(row['families']),

    output_line = [ str(row['chrom']),
                    str(row['start']),
                    str(row['end']),
                    str(row['ref']),
                    str(row['alt']),
                    str(row['rs_ids']),
                    #str(row['filter']),
                    str(row['cosmic_ids']),
                    #str(row['qual']),
                    #str(row['qual_depth']),
                    str(row['gene']),
                    str(row['transcript']),
                    str(row['num_hom_ref']),
                    str(row['num_het']),
                    str(row['num_hom_alt']),
                    #str(row['families']),
                    str(row['variant_samples']),
                    str(row['het_samples']),
                    str(row['hom_alt_samples']),
                    str(row['pfam_domain']),
                    str(row['codon_change']),
                    str(row['aa_change']),
                    str(row['impact']),
                    str(row['impact_severity']),
                    str(row['type']),
                    str(row['sub_type']),
                    str(row['call_rate']),
                    str(row['aaf']),
                    str(row['in_esp']),
                    str(row['aaf_esp_ea']),
                    str(row['aaf_esp_all']),
                    str(row['in_1kg']),
                    str(row['aaf_1kg_amr']),
                    str(row['aaf_1kg_eur']),
                    str(row['aaf_1kg_all']),
                    str(row['in_exac']),
                    str(row['aaf_exac_all']),
                    str(row['aaf_adj_exac_all']),
                    str(row['in_omim']),
                    str(row['clinvar_sig']),
                    str(row['clinvar_disease_name']),
                    str(row['clinvar_dbsource']),
                    str(row['clinvar_dbsource_id']),
                    str(row['clinvar_origin']),
                    str(row['clinvar_dsdb']),
                    str(row['clinvar_dsdbid']),
                    str(row['clinvar_disease_acc']),
                    str(row['clinvar_in_locus_spec_db']),
                    str(row['clinvar_on_diag_assay']),
                    str(row['gms_illumina']),
                    str(row['rmsk']),
                    str(row['in_cpg_island']),
                    str(row['in_segdup']),
                    str(row['is_conserved']),
                    str(row['in_cse']),
                    str(row['grc']),
                    str(row['gms_illumina']),
                    str(row['gerp_bp_score']),
                    str(row['gerp_element_pval']),
                    str(row['recomb_rate']),
                    str(row['cadd_raw']),
                    str(row['cadd_scaled']),
                    str(row['fitcons']) 
                  ]

    GType = row['gt_types']
    GT = row['gts']
    DP = row['gt_depths']
    GQ = row['gt_quals']
    R0 = row['gt_ref_depths']
    A0 = row['gt_alt_depths']

    for item in samples:
        #print item
        idx = smp2idx[item]
        #print idx
        if GT[idx] is None:
            sampleGT = "." 
        else:
            sampleGT = GT[idx]
            sampleDP = DP[idx]
            sampleGQ = GQ[idx]
            sampleR0 = R0[idx]
            sampleA0 = A0[idx]
            sampleGType = GType[idx]

        if ( float(sampleR0) + float(sampleA0) ) > 0 :
            sampleAF = float(sampleA0)/( float(sampleR0) + float(sampleA0) )
        else:
            sampleAF = 0

        metric = str(sampleGT) + ":" + str(sampleR0) + ":" + str(sampleA0) + ":" + str(sampleDP) + ":" + str(sampleGQ)
        #metric = str(sampleGT) + ":" + str(sampleR0) + ":" + str(sampleA0) + ":" + str(sampleDP)
        output_line.append(str(sampleDP))
        output_line.append(str(sampleAF))
        #output_line.append(str(sampleGType))
        output_line.append(str(metric))

    print "\t".join(output_line)

    #except KeyError:
    #    pass

	#output_line = [ str(row['chrom']),str(row['start']),str(row['end']),str(row['ref']),str(row['alt']),str(row['rs_ids']),str(row['cosmic_ids']),str(row['cosmic_ids']),str(row['qual']),str(row['qual_depth']),str(row['gene']),str(row['transcript']),str(row['num_hom_ref']),str(row['num_het']),str(row['num_hom_alt']),str(row['variant_samples']),str(row['HET_samples']),str(row['HOM_ALT_samples']),'None','None','None','None',str(row['pfam_domain']),str(row['codon_change']),str(row['aa_change']),str(row['impact']),str(row['impact_severity']),str(row['type']),str(row['sub_type']),str(row['call_rate']),str(row['aaf']),str(row['in_esp']),str(row['aaf_esp_ea']),str(row['aaf_esp_all']),str(row['in_1kg']),str(row['aaf_1kg_amr']),str(row['aaf_1kg_eur']),str(row['aaf_1kg_all']),str(row['in_exac']),str(row['aaf_exac_all']),str(row['aaf_adj_exac_all']),str(row['in_omim']),str(row['clinvar_sig']),str(row['clinvar_disease_name']),str(row['clinvar_dbsource']),str(row['clinvar_dbsource_id']),str(row['clinvar_origin']),str(row['clinvar_dsdb']),str(row['clinvar_dsdbid']),str(row['clinvar_disease_acc']),str(row['clinvar_in_locus_spec_db']),str(row['clinvar_on_diag_assay']),str(row['gms_illumina']),str(row['rmsk']),str(row['in_cpg_island']),str(row['in_segdup']),str(row['is_conserved']),str(row['in_cse']),str(row['cadd_raw']),str(row['cadd_scaled']),str(row['fitcons'])]
        
	#GT = row['gts']
	#DP = row['gt_depths']
    #GQ = row['gt_quals']
	#R0 = row['gt_ref_depths']
        #A0 = row['gt_alt_depths']

        #for item in samples:
        #       idx = smp2idx[item]
	    #	sampleGT = GT[idx]
        #       sampleDP = DP[idx]
        #       sampleGQ = GQ[idx]
	    #	sampleR0 = R0[idx]
        #       sampleA0 = A0[idx]

	    #	if ( float(sampleR0) + float(sampleA0) ) > 0 :
        #                sampleAF = float(sampleA0)/( float(sampleR0) + float(sampleA0) )
        #        else:
        #                sampleAF = 0

        #    metric = str(sampleGT) + ":" + str(sampleR0) + ":" + str(sampleA0) + ":" + str(sampleDP) + ":" + str(sampleGQ)
        #    output_line.append(str(sampleDP))
	    #	 output_line.append(str(sampleAF))
        #    output_line.append(str(metric))
        
        #print "\t".join(output_line)
