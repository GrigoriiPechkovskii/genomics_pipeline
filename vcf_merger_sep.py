import os
import pandas as pd

import vcf_merger
import pipeline_base


directory = os.getcwd()
file_vcf = directory + '/test/' + 'test_merger_alignment_checker.vcf'
file_fasta = directory + '/test/' +'test_merger.fna'

#log_file_path = directory + '/test/' + 'log.txt'
#log_file = open(log_file_path,'a')

fasta = pipeline_base.SequenceFasta(file_fasta)
fasta.seq_process(strip=True)
sequence = ''.join(fasta.seq_lst)


header,vcf_head,type_head = pipeline_base.vcf_head_process(file_vcf)


vcf = pd.read_csv(file_vcf,sep='\t',header = header,dtype=type_head)
vcf.index = vcf.index.astype(str) + '_' + vcf['#CHROM'].values.astype(str) + '_' +  vcf['POS'].values.astype(str)


vcf_merged = vcf_merger.merge_window([['A',2, 4]],vcf.copy(),ref_sequence=sequence,fullcheck=False,log_file=False)