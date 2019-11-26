print('start')
import os

import pandas as pd

file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker.vcf'
file_intervals1 = '/home/strain4/Desktop/genomics_pipline/interval1.txt'
file_intervals2 = '/home/strain4/Desktop/genomics_pipline/interval2.txt'

ref_assemble_name = 'A'

header = 5
vcf = pd.read_csv(file_vcf,sep='\t',header = header)#!15
vcf.index =vcf['POS'].values

#intervals_path
intervals_path = [i for i in os.listdir() if 'interval' in i]

def intervals_concat(intervals_path):
    '''intervals_path = list of path file interval'''
    intervals = []
    for interval in intervals_path:
        intervals.append(pd.read_csv(interval,sep='\t'))
    intervals_full = pd.concat(intervals,sort=False)

    return intervals_full

def vcf_correcter(vcf,intervals_full):

    assemble_names = list(vcf.columns[9:])
    print(assemble_names)
    for pos in vcf.index:
        pos_end = pos + len(vcf.loc[pos,'REF']) - 1
        #intervals_slice_contain = intervals_full[(intervals_full[ref_assemble_name + '_start'] <= pos) & (pos_end <= intervals_full[ref_assemble_name + '_end'])]
        intervals_slice_contain = intervals_full[( pos_end >=  intervals_full[ref_assemble_name + '_start']) & (pos <= intervals_full[ref_assemble_name + '_end'])]
        for assemble_name in assemble_names:
            print(pos,pos_end,assemble_name)
            #print(intervals_slice_contain)
            if vcf.loc[pos,assemble_name] == '.':
                print(intervals_slice_contain)
                #print(intervals_slice_contain[assemble_name + '_start'].notna(),intervals_slice_contain[assemble_name + '_start'] != 0)
                print(all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0))

                if assemble_name + '_start' not in intervals_slice_contain.columns:
                    #assert assemble_name + '_start' in intervals_slice_contain.columns, 'Warning assemble ' + assemble_name  + ' not in intervals'
                    continue
                    print('Warning assemble',assemble_name,'not in intervals')
                #elif  any(intervals_slice_contain[assemble_name + '_start'].notna()) or any(intervals_slice_contain[assemble_name + '_start'] != 0):#!!!!!
                elif all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0):
                    vcf.at[pos,assemble_name] = '0'
    return vcf
        
intervals_full = intervals_concat(intervals_path)

vcf = vcf_correcter(vcf,intervals_full)

vcf.to_csv('/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker2.vcf',sep='\t',index=False)

print('end')