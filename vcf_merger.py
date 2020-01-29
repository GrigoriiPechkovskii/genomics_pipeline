#!/usr/bin/python3

#By Grigorii Pechkovskii

import os
import re
import argparse
import time
import timeit

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', 25)

print('start')


parser = argparse.ArgumentParser()

parser.add_argument('-v', '--vcf',action='store', help='File vcf')
parser.add_argument('-r', '--ref',action='store', help='Reference fasta')
parser.add_argument('-g', '--gbk-file',action='store', help='File gbk')
parser.add_argument('-d', '--dir',action='store', help='Work directory')
parser.add_argument('-l', '--log',action='store', help='Log file')
parser.add_argument('-i', '--interval',action='store', help='Interval path')
parser.add_argument('-t', '--header-vcf',action='store',default='15', help='Header vcf')
parser.add_argument('-o', '--out-file',action='store',default=os.getcwd()+'/test_merged.vcf' , help='Out directory')

#parser.add_argument('-c', '--contig',action='store', help='Contig')
#parser.add_argument('-n', '--name-vcf',action='store',default='test_merged.vcf', help='Name vcf')

file_vcf = parser.parse_args().vcf
file_fasta = parser.parse_args().ref
file_gbk = parser.parse_args().gbk_file
work_dir = parser.parse_args().dir
log_file_path = parser.parse_args().log
header = int(parser.parse_args().header_vcf)
out_file = parser.parse_args().out_file

pd.set_option('display.max_columns', 10)


local = False
test = True

if not(local or test):
    intervals_path = parser.parse_args().interval.split(',')
    print('intervals_path=',intervals_path)
    ref_assemble_name = os.path.basename(file_fasta)[:-4]

if local:
    directory = os.getcwd()
    file_vcf = '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/merged.vcf'
    file_fasta = '/home/strain4/Desktop/fin_script/test_genomics_pipline/genome/GCF_000008445.1_ASM844v1_genomic.fna'
    file_gbk = directory + '/test/' + 'AmesAncestor_GCF_000008445.1.gbk'
    work_dir = '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out'
    out_file = '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/vcf_merged_5555.vcf'
    log_file_path = work_dir + '/' + 'log_test4444.txt'
    header = 15

    #intervals_path = [work_dir+'/'+i for i in os.listdir(work_dir) if 'interval' in i]
    intervals_path = ['/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_0.bed', '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_1.bed']

    ref_assemble_name = os.path.basename(file_fasta)[:-4]

if test:
    directory = os.getcwd()
    file_fasta = directory + '/test/' +'test_merger.fna'
    file_gbk = directory + '/test/' + 'AmesAncestor_GCF_000008445.1.gbk'
    file_vcf = directory + '/test/' + 'test_merger_alignment_checker.vcf'
    work_dir = directory + '/test'
    log_file_path = directory + '/test/' + 'log.txt'
    header = 5
    ref_assemble_name = 'A'
    out_file = work_dir + '/' + 'test_merger.vcf'
    intervals_path = [work_dir+'/'+i for i in os.listdir(work_dir) if '.bed' in i]
    #intervals_path = [work_dir+'/'+i for i in os.listdir(work_dir) if 'interval' in i]

log_file = open(log_file_path,'a')

vcf = pd.read_csv(file_vcf,sep='\t',header = header)#!15
vcf.index =vcf['POS'].values

class SequenceFasta():
    '''Processes a fasta file'''
    def __init__(self,file_dir):
        self.file_dir=file_dir
        self.name_lst = []
        self.seq_lst = []
    def just_named(self):
        '''Works only with sequence names'''
        with open(self.file_dir) as file_opened:
            for line in file_opened:
                if '>' in line:
                    print(line)
                    self.name_lst.append(line)
    def seq_process(self,strip=False):
        seq_tmp = ''
        with open(self.file_dir) as file_opened:            
            for line in file_opened:
                if '>' in line:
                    self.seq_lst.append(seq_tmp)
                    seq_tmp = ''
                    if strip:
                        self.name_lst.append(line.rstrip())
                    else:
                        self.name_lst.append(line)  
                else:
                    if strip:
                        seq_tmp += line.rstrip()
                    else:
                        seq_tmp += line
        self.seq_lst.append(seq_tmp)
        del self.seq_lst[0]
        self.seq_len = len(self.seq_lst)

fasta = SequenceFasta(file_fasta)
fasta.seq_process(strip=True)
sequence = ''.join(fasta.seq_lst)

def merge_window(intervals,vcf,ref_sequence=sequence,fullcheck=True,ignored=True):
    '''Function that joins a given window in a vcf file by position'''
    global vcf_slice
    intervals_used = []
    for pos_start,pos_end in intervals:

        flag_ignored = False        

        if any([pos_start_used <= pos_start <= pos_end_used for pos_start_used,pos_end_used in intervals_used]):
            print('Warning this interval used, interval ignored',[pos_start,pos_end])
            log_file.write('Warning this interval used, interval ignored ' + str([pos_start,pos_end])+'\n')
            continue
                
        contig = vcf.loc[pos_start,'#CHROM']

        vcf_slice = vcf[(pos_start <= vcf['POS']) & (vcf['POS'] <= pos_end) & (vcf['#CHROM']==contig)]       

        if vcf_slice.shape[0] <= 1:
            continue

        sample_bin_dict = {i:None for i in vcf_slice.iloc[:,9:]}
        position_first = int(vcf_slice.iloc[0]['POS'])
        position_start = int(vcf_slice.iloc[0]['POS'])
        max_over = max([len(ref) + pos - 1 for ref,pos in zip(vcf_slice['REF'],vcf_slice.index)])

        sample_dict = {i:list(sequence[position_first-1:max_over]) for i in vcf_slice.iloc[:,9:]}
        sample_dict_mem = {i:position_first for i in vcf_slice.iloc[:,9:]}

        ref_mod = sequence[position_first-1:max_over]
        sample_will_deleted = set()
        
        #check for if overlay out interval
        if fullcheck:
            '''
            vcf_slice_end = vcf[(vcf['POS'] > pos_end) & (vcf['#CHROM']==contig)]
            
            if not vcf_slice_end.empty:  
                position_outside = int(vcf_slice_end.iloc[0]['POS'])
                for position_row in vcf_slice.index[:-1]:
                    if position_row + len(vcf_slice.loc[position_row]['REF']) - position_outside > 0:
                        print('Warning overlay window','interval',[pos_start,pos_end])
                        log_file.write('Warning overlay window, interval '+str([pos_start,pos_end])+'\n')
            '''
            flag = False
            if  any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                print('Warning indeterminate genotype of the sample1 interval ignored',[pos_start,pos_end])
                log_file.write('Warning indeterminate genotype of the sample1 interval ignored ' + str([pos_start,pos_end])+'\n')

            for sample in sample_dict:
                if  not all(vcf_slice[sample].isin(['.'])) and not all(~vcf_slice[sample].isin(['.'])):#!
                    #print('Warning indeterminate genotype of the sample')
                    #assert False, "'Warning indeterminate genotype of the sample'"
                    flag = True
            if flag:
                #print('Warning2 indeterminate genotype of the sample, interval ignored',[pos_start,pos_end])
                #log_file.write('Warning2 indeterminate genotype of the sample, interval ignored ' + str([pos_start,pos_end])+'\n')
                if ignored:
                    print('Warning2 indeterminate genotype of the sample, interval ignored',[pos_start,pos_end])
                    log_file.write('Warning2 indeterminate genotype of the sample, interval ignored ' + str([pos_start,pos_end])+'\n')
                    continue

            if not (vcf_slice.iloc[:,9:]).isin(['.']).any().any():#!req to change
                intervals_slice_contain = intervals_alignment[( pos_end >=  intervals_alignment[ref_assemble_name + '_start']) & (pos_start <= intervals_alignment[ref_assemble_name + '_end'])]
                if intervals_slice_contain.isin([0]).any().any(): 
                    if ignored:
                        print('Warning gap between defined variant, interval ignored',[pos_start,pos_end])
                        log_file.write('Warning gap between defined variant, interval ignored ' + str([pos_start,pos_end])+'\n')
                        continue
        variant = []        
        for position_row in vcf_slice.index:
            variant = [vcf_slice.loc[position_row]['REF']] + vcf_slice.loc[position_row]['ALT'].split(',')
            position_inter = int(vcf_slice.loc[position_row]['POS'])
            #vcf.loc[position_row,'REF'] = ref_mod #maybe need
            for sample in sample_dict:
                if vcf_slice.loc[position_row][sample] == '.':
                    sample_dict[sample]
                    sample_will_deleted.add(sample)
                else:
                    if int(vcf_slice.loc[position_row][sample]) != 0:
                        position_row_shift = position_row - position_first
                        var_len = len(variant[int(vcf_slice.loc[position_row][sample])])
                        ref_len = len(vcf_slice.loc[position_row]['REF'])
                        
                        if any([i=='' or type(i)==list for i in sample_dict[sample][position_row_shift:position_row_shift+ref_len]]):
                            print('Error, position was used for variant, maybe wrong output,interval ignored',[pos_start,pos_end])
                            log_file.write('Error, position was used for variant, maybe wrong output,interval ignored ' + str([pos_start,pos_end])+'\n')
                            flag_ignored = True
                            break

                        #for cut end sequence maybe
                        ref_var = vcf_slice.loc[position_row]['REF']
                        alt_var = variant[int(vcf_slice.loc[position_row][sample])]
                        #if len(ref_var)>1 or len(alt_var)>1:
                        alt_cut = 0
                        '''
                        for sym_ref, sym_alt in zip(ref_var[::-1],alt_var[::-1]):
                            if sym_ref == sym_alt:
                                pass
                                #ref_len -= 1
                                #alt_cut -= 1
                                #alt_var = alt_var[:-1]
                            else:
                                break
                        '''
                        alt_cut = None if alt_cut==0 else alt_cut


                        if len(variant[int(vcf_slice.loc[position_row][sample])]) == ref_len:#!
                            sample_dict[sample][position_row_shift:position_row_shift+ref_len] =  ([list(variant[int(vcf_slice.loc[position_row][sample])])])[:alt_cut]

                        elif len(variant[int(vcf_slice.loc[position_row][sample])]) != ref_len:
                            sample_dict[sample][position_row_shift:position_row_shift+ref_len] = ([list(variant[int(vcf_slice.loc[position_row][sample])])])[:alt_cut] + [''] * (ref_len - 1) #(len(sample_dict[sample][position_row_shift:position_row_shift+ref_len]) - 1)
                            
                        sample_dict_mem[sample] += len(vcf_slice.loc[position_row]['REF'])
                           
                    else:
                        continue
            position_start = position_inter       
        
        if flag_ignored:
            continue
        intervals_used.append([pos_start,pos_end])

        #for redefine bin in sample
        bin_value = 1
        used = []
        used_bin = []
        flag = False
        sample_lst = []

        for sample in sample_will_deleted:
            del sample_dict[sample]
            sample_bin_dict[sample] = '.'


        for sample in sample_dict:
            sample_dict[sample] = [i if type(i) == str else ''.join(i) for i in sample_dict[sample]]
        sample_dict = {i:''.join(j) for i,j in sample_dict.items()}

        for sample in sample_dict:
            if all(vcf_slice[sample].astype(int) == 0): #!need check . values        
                sample_bin_dict[sample] = 0
            else:
                if sample not in used:
                    for sample_compared in sample_dict:
                        if (vcf_slice.iloc[:,9:][sample_compared].astype(int)).equals(vcf_slice.iloc[:,9:][sample].astype(int)) and sample_compared not in used:# and list(vcf_slice.iloc[:,9:][sample_compared].values) not in used_bin:# and sample_compared != sample:
                            sample_bin_dict[sample_compared] = bin_value
                            used += [sample_compared]
                            used_bin.append(list(vcf_slice.iloc[:,9:][sample_compared].values))
                            flag = True
                if flag == True:

                    sample_lst += [sample_dict[sample]]
                    bin_value += 1
                    flag = False

        alt = ','.join(sample_lst)
        info = 'INDEL'
        print('Add compressed interval to vcf',pos_start,pos_end)
        log_file.write('Add compressed interval to vcf ' + str([pos_start,pos_end])+'\n')

        vcf_variant_dict = {'#CHROM':contig,'POS':position_first,'ID':'.','REF':ref_mod,
                    'ALT':alt,'QUAL':40,'FILTER':'PASS','INFO':info,'FORMAT':'GT'}
        vcf_variant_dict.update(sample_bin_dict)

        for_merged = pd.Series(vcf_variant_dict,name=position_first)
        vcf.drop(vcf_slice.index,inplace=True)
        vcf = vcf.append(for_merged,ignore_index=False)#!True

    vcf.sort_values(by=['#CHROM','POS'],ascending=[False,True],inplace=True)

    return vcf

def definer_overlap_window(vcf):
    over = 0
    interval_noexact = []
    interval_exact = []
    pos = vcf.iloc[0]["POS"]
    for num,position_row in enumerate(vcf.loc[pos:].index):
        contig = vcf.loc[position_row,'#CHROM']
        interval = []
        interval2 = []
        window_sum = position_row + len(vcf.loc[position_row]['REF'])-1# + 5
        #print(window_sum)
        for position_row2 in vcf.index[num+1:]:
            if position_row2 <= window_sum:
                over = position_row2+len(vcf.loc[position_row2]['REF'])-1# + 5
                #print('minus=',window_sum,over,window_sum - over)
                if over > window_sum: #or 
                    window_sum = over 

                if window_sum - position_row>1000:
                    print('Warning interval contain large indel','interval ignored',[position_row,position_row2],'Indel lenght =',window_sum-position_row)#for filter large variant
                    log_file.write('Warning interval contain large indel, interval ignored ' + str([position_row,position_row2])+' Indel lenght = ' + str(window_sum-position_row) +'\n')
                    break

                #vcf_slice = vcf[(position_row <= vcf['POS']) & (vcf['POS'] <= position_row2) & (vcf['#CHROM']==contig)]#!window_sum
                #if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                    #print('Warning interval contain indeterminate variant,interval ignored',[position_row,position_row2])
                    #log_file.write('Warning interval contain indeterminate variant,interval ignored ' + str([position_row,position_row2])+'\n')
                    #break
                interval = [position_row,window_sum]
                interval2 = [position_row,position_row2]
            else:
                break
        
        if interval:
            interval_noexact.append(interval)
            
        if interval2:
            interval_exact.append(interval2)

    return interval_exact,interval_noexact

def intervals_concat(intervals_path):
    '''intervals_path = list of path file interval'''
    intervals = []
    for interval in intervals_path:

        if test:
            intervals.append(pd.read_csv(interval,sep='\t').fillna(0))#!!index col        
        else:
            intervals.append(pd.read_csv(interval,sep='\t',index_col=0).fillna(0))#!!index col

    intervals_full = pd.concat(intervals,sort=False)

    return intervals_full


def intervals_concat_bed(intervals_path):
    '''intervals_path = list of path file interval'''
    intervals = pd.DataFrame()
    intervals_full = intervals.append([pd.read_csv(path_bed,sep='\t') for path_bed in intervals_path],ignore_index=True,sort=False)

    return intervals_full

time_block_start = time.time()
time_initial = time.time()
def timecheck(name = 'name_check'):
    '''TimeCheck for testing'''
    global time_block_start
    check_block_sec = round(time.time() - time_block_start,3)
    check_block_min = round(check_block_sec / 60,3)

    check_initial_sec = round(time.time() - time_initial,3)
    check_initial_min = round(check_initial_sec / 60,3)       
    report_line = '{0:.<30}{1:.<15}{2:.<10}{3:.<10}{4:.<10}'.format(name,check_block_sec,check_block_min,check_initial_sec,check_initial_min) + '\n'
    print(report_line,end='')
    #with open(HOME_DIRECTORY + 'time_report.txt','a') as time_report:
    #    time_report.write(report_line)    
    time_block_start = time.time()
    return 

def vcf_corrector_1ver(vcf,intervals_alignment):
    sum_0 = 0
    sum_NA = 0
    assemble_names = list(vcf.columns[9:])
    for assemble_name in assemble_names:
        assemble_slice = vcf[vcf[assemble_name].astype(str).eq('.')]
        assemble_slice_full  = vcf[assemble_name].copy()
        assemble_slice_pos = assemble_slice[assemble_name].copy()
        pos_end_array = assemble_slice[assemble_name].index + assemble_slice['REF'].map(len) - 1
        for pos in assemble_slice_pos.index:
            pos_end = pos_end_array[pos]
            intervals_slice_contain = intervals_alignment[( pos_end >=  intervals_alignment[ref_assemble_name + '_start']) & (pos <= intervals_alignment[ref_assemble_name + '_end'])]
            if assemble_name + '_start' not in intervals_slice_contain.columns:
                #assert assemble_name + '_start' in intervals_slice_contain.columns, 'Warning assemble ' + assemble_name  + ' not in intervals'
                #print('Error assemble',assemble_name,'not in intervals')
                log_file.write('Error assemble '+assemble_name+' not in intervals'+'\n')
                continue
            if all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0):
                assemble_slice_full[pos] = 0
                sum_0 += 1
            else:
                sum_NA += 1

        vcf[assemble_name] = assemble_slice_full
    return vcf

def vcf_corrector_2ver(vcf,intervals_alignment):
    sum_0 = 0
    sum_NA = 0
    assemble_names = list(vcf.columns[9:])
    for assemble_name in assemble_names:
        assemble_slice = vcf[assemble_name][vcf[assemble_name].astype(str).eq('.')]
        assemble_slice_full  = vcf[assemble_name].copy()
        pos_end_array = assemble_slice.index + vcf.loc[assemble_slice.index]['REF'].map(len) - 1
        for pos in assemble_slice.index:
            pos_end = pos_end_array[pos]
            intervals_slice_contain = intervals_alignment[( pos_end >=  intervals_alignment[ref_assemble_name + '_start']) & (pos <= intervals_alignment[ref_assemble_name + '_end'])]
            if assemble_name + '_start' not in intervals_slice_contain.columns:
                #assert assemble_name + '_start' in intervals_slice_contain.columns, 'Warning assemble ' + assemble_name  + ' not in intervals'
                #print('Error assemble',assemble_name,'not in intervals')
                log_file.write('Error assemble '+assemble_name+' not in intervals'+'\n')
                continue
            if all(intervals_slice_contain[assemble_name + '_start'][intervals_slice_contain[assemble_name + '_start'].notna()] != 0):
                assemble_slice_full[pos] = 0
                sum_0 += 1
            else:
                sum_NA += 1
        vcf[assemble_name] = assemble_slice_full
    return vcf


def vcf_corrector_bed(vcf,intervals_alignment):
    sum_0 = 0
    sum_NA = 0
    assemble_names = list(vcf.columns[9:])
    for index,bedline in intervals_alignment.iterrows():
        if bedline['name'] in vcf.columns:
            if (bedline['start_position_ref']!=0 or bedline['end_position_ref']!=0) and (bedline['start_position_alt']!=0 or bedline['end_position_alt']!=0):
                vcf.loc[(vcf['POS']>=bedline['start_position_ref']) & (vcf['POS']<=bedline['end_position_ref']) & (vcf[bedline['name']]=='.'),bedline['name']] = 0
                sum_0+=1
            else:
                sum_NA+=1
        else:
            print('Error sample not in vcf')
            
    return vcf


def contig_finder_gbk(file_gbk_dir):
    ''' '''
    with open(file_gbk) as file_gbk_opened:
        file_gbk_read = file_gbk_opened.read()
        find_locus = re.findall(r'LOCUS\s+(.*?)\s\s+',file_gbk_read)
        find_source = re.findall(r'\s\s+source\s+(.*?)\s\s+',file_gbk_read)

    find_source = [[*map(int,(i.split('..')))] for i in find_source]
    find_source_real = find_source.copy()
    for source_num in range(1,len(find_source)):
        find_source[source_num] = [find_source[source_num][0] + find_source[source_num-1][1],
                                   find_source[source_num][1] + find_source[source_num-1][1]]
    return find_locus, find_source,find_source_real

def contig_definder(position,find_locus,find_source): 
    ''' ''' 
    for locus,source in zip(find_locus,find_source):
        if (source[0] <= position <= source[1]):
            position_real = position - source[0]+ 1#!
            return locus,position_real

def position_editer(vcf):
    pos_lst = []
    for position in vcf['POS']:
        contig,position_real = contig_definder(position,find_locus,find_source)
        pos_lst.append(position_real)
        #vcf.loc[position,'POS'] = position_real
    vcf['POS'] = np.array(pos_lst)
    vcf.index = vcf['POS']
    return vcf

vcf = pd.read_csv(file_vcf,sep='\t',header = header)
vcf.index =vcf['POS'].values

log_file = open(log_file_path,'a')
log_file.write('\n'+'Start vcf_corrector'+'\n\n')
print('\n'+'Start intervals_concat'+'\n\n')
#intervals_alignment = intervals_concat(intervals_path)
intervals_alignment = intervals_concat_bed(intervals_path)

timecheck('first')
print('\n'+'Start vcf_corrector'+'\n\n')
#vcf_correct_bed = vcf_corrector_bed(vcf.copy(),intervals_alignment)
timecheck('vcf_corrector_bed_before')
for i in range(1):
    vcf_correct_bed = vcf_corrector_bed(vcf.copy(),intervals_alignment)
timecheck('vcf_corrector_bed_after')
'''
intervals_alignment = intervals_concat(intervals_path)
timecheck('vcf_corrector1_before')
for i in range(1):
    vcf_correct = vcf_corrector_1ver(vcf.copy(),intervals_alignment)
timecheck('vcf_corrector1_after')

timecheck('vcf_corrector2_before')
for i in range(1):
    vcf_correct = vcf_corrector_2ver(vcf.copy(),intervals_alignment)
timecheck('vcf_corrector2_after')
'''

log_file = open(log_file_path,'a')
log_file.write('\n'+'Start definer_overlap_window'+'\n\n')
print('\n'+'Start definer_overlap_window'+'\n\n')
interval_exact,interval_noexact = definer_overlap_window(vcf_correct_bed)
#print('interval_noexact =',interval_noexact,'interval_exact=',interval_exact)

log_file = open(log_file_path,'a')
log_file.write('\n'+'Start merge_window'+'\n\n')
print('\n'+'Start merge_window'+'\n\n')
vcf_merged = merge_window(interval_exact,vcf_correct_bed.copy(),fullcheck=False)
#vcf_merged = merge_window([[226253, 226261],[85540,85551]],vcf_correct.copy())

find_locus, find_source ,find_source_real = contig_finder_gbk(file_gbk)
vcf_merged = position_editer(vcf_merged.copy())

vcf_merged.to_csv(out_file,sep='\t',index=False)

log_file.close()
timecheck('end')
print('end')