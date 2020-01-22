print('start')
import os

import pandas as pd

file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker.vcf'
file_intervals1 = '/home/strain4/Desktop/genomics_pipline/interval1.txt'
file_intervals2 = '/home/strain4/Desktop/genomics_pipline/interval2.txt'
work_dir = '.'
ref_assemble_name = 'A'
header = 5


file_vcf = '/home/strain4/Desktop/test_genomics_pipline/exp4_test/vcf_out/merged.vcf'
work_dir = '/home/strain4/Desktop/test_genomics_pipline/exp4_test/vcf_out'
ref_assemble_name = 'GCF_000008445.1_ASM844v1_genomic'
header = 15

vcf = pd.read_csv(file_vcf,sep='\t',header = header)#!15
vcf.index =vcf['POS'].values


#intervals_path
intervals_path = [work_dir+'/'+i for i in os.listdir(work_dir) if 'interval' in i]

def intervals_concat(intervals_path):
    '''intervals_path = list of path file interval'''
    intervals = []
    for interval in intervals_path:
        intervals.append(pd.read_csv(interval,sep='\t'))
    intervals_full = pd.concat(intervals,sort=False)

    return intervals_full


def vcf_correcter(vcf,intervals_full):
    sum_0 = 0
    sum_NA = 0
    assemble_names = list(vcf.columns[9:])
    #print(assemble_names)
    for pos in vcf.index:
        pos_end = pos + len(vcf.loc[pos,'REF']) - 1
        #intervals_slice_contain = intervals_full[(intervals_full[ref_assemble_name + '_start'] <= pos) & (pos_end <= intervals_full[ref_assemble_name + '_end'])]
        intervals_slice_contain = intervals_full[( pos_end >=  intervals_full[ref_assemble_name + '_start']) & (pos <= intervals_full[ref_assemble_name + '_end'])]
        for assemble_name in assemble_names:
            #print(pos,pos_end,assemble_name)
            #print(intervals_slice_contain)
            if vcf.loc[pos,assemble_name] == '.':
                #print(intervals_slice_contain)
                #print(intervals_slice_contain[assemble_name + '_start'].notna(),intervals_slice_contain[assemble_name + '_start'] != 0)
                #print(all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0))

                if assemble_name + '_start' not in intervals_slice_contain.columns:
                    #assert assemble_name + '_start' in intervals_slice_contain.columns, 'Warning assemble ' + assemble_name  + ' not in intervals'
                    print('Warning assemble',assemble_name,'not in intervals')
                    continue
                #elif  any(intervals_slice_contain[assemble_name + '_start'].notna()) or any(intervals_slice_contain[assemble_name + '_start'] != 0):#!!!!!
                if all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0):
                    vcf.at[pos,assemble_name] = '0'
                    sum_0 += 1
                else:
                    sum_NA += 1
    print('sum_0  = ',sum_0,'\n','sum_NA = ',sum_NA,sep='')
    return vcf
        
intervals_full = intervals_concat(intervals_path)

vcf = vcf_correcter(vcf,intervals_full)

vcf.to_csv('/home/strain4/Desktop/test_genomics_pipline/test_merger.vcf',sep='\t',index=False)

print('end')