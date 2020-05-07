#!/usr/bin/python3

#By Grigorii Pechkovskii

import os
import re
import argparse
import time
import timeit

import numpy as np
import pandas as pd

import pipeline_base


def merge_window(intervals,vcf,ref_sequence,log_file=False,fullcheck=True,ignored=True,info_just_indel=False,drop_info=False,intervals_alignment_bool=False,intervals_alignment=None):
    '''Function that joins a given window in a vcf file by position'''
    global vcf_slice
    intervals_used = []
    for contig,pos_start,pos_end in intervals:

        flag_ignored = False        

        if any([(pos_start_used <= pos_start <= pos_end_used) and contig==contig_used for contig_used,pos_start_used,pos_end_used in intervals_used]):
            print('Warning this interval used, interval ignored',[contig,pos_start,pos_end])
            
            if log_file: 
                log_file.write('Warning this interval used, interval ignored ' + str([contig,pos_start,pos_end])+'\n')
            
            continue
                
        #contig = vcf.loc[pos_start,'#CHROM']
        #print(contig,pos_start,pos_end)
        vcf_slice = vcf[(pos_start <= vcf['POS']) & (vcf['POS'] <= pos_end) & (vcf['#CHROM']==contig)]

        #print(vcf_slice)       

        if vcf_slice.shape[0] <= 1:
            print('Warning vcf_slice <= 1',contig,pos_start,pos_end)

            if log_file: 
                log_file.write('Warning vcf_slice <= 1' +'\n')

            continue

        if drop_info:
            if len(set(vcf_slice['INFO'].values) - set(drop_info)) == 0:
                print('Warning vcf slice contain only', *drop_info,set(vcf_slice['INFO'].values),contig,pos_start,pos_end)
                #print('Warning vcf slice contain only',str([contig,pos_start,pos_end]))

                if log_file:
                    log_file.write('Warning vcf slice contain only' + str([contig,pos_start,pos_end])+'\n')

                continue


        #sample_bin_dict = {i:None for i in vcf_slice.iloc[:,9:]}#!!!maybe need
        sample_bin_dict = dict()

        position_first = int(vcf_slice.iloc[0]['POS'])
        index_first = vcf_slice.index[0]
        position_start = int(vcf_slice.iloc[0]['POS'])
        max_over = max([len(ref) + pos - 1 for ref,pos in zip(vcf_slice['REF'],vcf_slice['POS'])])

        sample_dict = {i:list(ref_sequence[position_first-1:max_over]) for i in vcf_slice.iloc[:,9:]}
        sample_dict_info = {i:'' for i in vcf_slice.iloc[:,9:]}

        sample_dict_mem = {i:position_first for i in vcf_slice.iloc[:,9:]}

        ref_mod = ref_sequence[position_first-1:max_over]
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
                
                if log_file:
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
                    
                    if log_file:
                        log_file.write('Warning2 indeterminate genotype of the sample, interval ignored ' + str([pos_start,pos_end])+'\n')
                    
                    continue

            if not (vcf_slice.iloc[:,9:]).isin(['.']).any().any():#!req to change
                intervals_slice_contain = intervals_alignment[( pos_end >=  intervals_alignment[ref_assemble_name + '_start']) & (pos_start <= intervals_alignment[ref_assemble_name + '_end'])]
                if intervals_slice_contain.isin([0]).any().any(): 
                    if ignored:
                        print('Warning gap between defined variant, interval ignored',[pos_start,pos_end])
                        
                        if log_file:
                            log_file.write('Warning gap between defined variant, interval ignored ' + str([pos_start,pos_end])+'\n')
                        
                        continue

        

        if intervals_alignment_bool:
            intervals_alignment_slice = intervals_alignment[(intervals_alignment['#contig'] == contig ) &
                #(intervals_alignment['name'] == sample ) &
                ((intervals_alignment['start_position_ref'] <= pos_start) &
                (intervals_alignment['end_position_ref']  >= pos_start))|

                ((intervals_alignment['start_position_ref'] <= pos_end)&
                   (intervals_alignment['end_position_ref'] >= pos_end))]
            
            sample_tmp_storage = []
            #print(pos_start,pos_end)
            for sample in sample_dict:
                #intervals_alignment['#contig'] == contig

                #print(sample,pos_start,pos_end)
                #print(intervals_alignment_slice[intervals_alignment_slice['name'] == sample])

                if intervals_alignment_slice[intervals_alignment_slice['name'] == sample].shape[0] >= 2:
                    sample_will_deleted.add(sample)
                    sample_tmp_storage.append(sample)
                    #print(intervals_alignment_slice[intervals_alignment_slice['name'] == sample])
                    #print('\n\n\n')           

            if sample_tmp_storage:
            	print('Warning merge different contig',contig,pos_start,pos_end,sample_tmp_storage)
            	if log_file:
            		log_file.write('Warning merge different contig ' + str([contig,pos_start,pos_end,sample_tmp_storage])+'\n')


        
        variant = []
        flag_break = False        
        for row_index in vcf_slice.index:
            if flag_break:
                break
            #position_row
            #print(row_index)
            variant = [vcf_slice.loc[row_index]['REF']] + vcf_slice.loc[row_index]['ALT'].split(',')
            position_inter = int(vcf_slice.loc[row_index]['POS'])
            #vcf.loc[position_row,'REF'] = ref_mod #maybe need
            for sample in sample_dict:


                if vcf_slice.loc[row_index][sample] == '.':
                    #sample_dict[sample]
                    sample_will_deleted.add(sample)
                else:
                    if int(vcf_slice.loc[row_index][sample]) != 0:
                        position_row_shift = int(vcf_slice.loc[row_index]['POS']) - position_first
                        var_len = len(variant[int(vcf_slice.loc[row_index][sample])])
                        ref_len = len(vcf_slice.loc[row_index]['REF'])

                        #print(sample)
                        #print(sample_dict,file=log_file)
                        
                        if any([i=='' or type(i)==list for i in sample_dict[sample][position_row_shift:position_row_shift+ref_len]]):
                            print(row_index,sample,'Error, position was used for variant, maybe wrong output,interval ignored',[pos_start,pos_end])
                            print('\n',position_row_shift,position_row_shift+ref_len,sample_dict,file=log_file)
                            
                            if log_file:
                                log_file.write('Error, position was used for variant, maybe wrong output,interval ignored ' + str([pos_start,pos_end])+'\n')
                            
                            flag_ignored = True
                            flag_break=True
                            break

                        #for cut end ref_sequence maybe
                        ref_var = vcf_slice.loc[row_index]['REF']
                        alt_var = variant[int(vcf_slice.loc[row_index][sample])]
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


                        if len(variant[int(vcf_slice.loc[row_index][sample])]) == ref_len:#!
                            sample_dict[sample][position_row_shift:position_row_shift+ref_len] =  ([list(variant[int(vcf_slice.loc[row_index][sample])])])[:alt_cut]
                            #sample_dict_info[sample] = sample_dict_info[sample]+vcf_slice.loc[row_index]['INFO'] + '_'

                        elif len(variant[int(vcf_slice.loc[row_index][sample])]) != ref_len:
                            sample_dict[sample][position_row_shift:position_row_shift+ref_len] = ([list(variant[int(vcf_slice.loc[row_index][sample])])])[:alt_cut] + [''] * (ref_len - 1) #(len(sample_dict[sample][position_row_shift:position_row_shift+ref_len]) - 1)
                        
                        sample_dict_info[sample] = sample_dict_info[sample]+vcf_slice.loc[row_index]['INFO'] + '+'

                        sample_dict_mem[sample] += len(vcf_slice.loc[row_index]['REF'])
                           
                    else:
                        continue
            position_start = position_inter       
        
        if flag_ignored:
            continue
        intervals_used.append([contig,pos_start,pos_end])

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
        

        sample_dict_info = {i:j[:-1] for i,j in sample_dict_info.items()}
        ref_mod
        #sample_bin_dict = dict()
        sample_seq_set = set()
        sample_uniq_dict = dict()
        info_lst = []
        sample_uniq_dict[ref_mod] = 0
        sample_seq_set.add(ref_mod)

        

        bin_value = 1
        for sample,sample_seq_uniq in sample_dict.items():
            if sample_seq_uniq not in sample_seq_set:
                sample_seq_set.add(sample_seq_uniq)
                sample_lst.append(sample_seq_uniq)

                info_lst.append(sample_dict_info[sample])

                sample_uniq_dict[sample_seq_uniq] = bin_value
                bin_value+=1

        for sample in sample_dict:
            sample_bin_dict[sample] = sample_uniq_dict[sample_dict[sample]]
        info = '&'.join(info_lst)

        #print('sample_bin_dict =',sample_dict)

        '''
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
        '''

        alt = ','.join(sample_lst)

        #print(vcf_slice['INFO'].values)
        
        if info_just_indel:
            info = 'INDEL'
        print('Add compressed interval to vcf',pos_start,pos_end,str(pos_end-pos_start))
        
        if log_file:
            log_file.write('Add compressed interval to vcf ' + str([pos_start,pos_end])+ str(pos_end-pos_start) + '\n')

        vcf_variant_dict = {'#CHROM':contig,'POS':position_first,'ID':'.','REF':ref_mod,
                    'ALT':alt,'QUAL':40,'FILTER':'PASS','INFO':info,'FORMAT':'GT'}
        vcf_variant_dict.update(sample_bin_dict)

        for_merged = pd.Series(vcf_variant_dict,name=index_first)
        vcf.drop(vcf_slice.index,inplace=True)
        vcf = vcf.append(for_merged,ignore_index=False)#!True

    vcf.sort_values(by=['#CHROM','POS'],ascending=[False,True],inplace=True)

    return vcf

def definer_overlap_window(vcf,overlap_extra=0,type_merge='hard',log_file=False):
    over = 0
    interval_noexact = []
    interval_exact = []
    pos = vcf.iloc[0]["POS"]
    #for num,position_row in enumerate(vcf.loc[pos:].index):
    for num,position_row in enumerate(vcf.index):
        contig = vcf.iloc[num]['#CHROM']

        interval = []
        interval_dict = {}
        interval2 = []
        window_sum = vcf.loc[position_row]['POS'] + len(vcf.loc[position_row]['REF'])-1 + overlap_extra
        #print(position_row,window_sum)
        for position_row2 in vcf.index[num+1:]:
            if vcf.loc[position_row2]['POS'] <= window_sum:
                
                if vcf.loc[position_row2]['#CHROM'] != contig:
                    print('Warning different contig',vcf.loc[position_row2]['#CHROM'], contig,position_row,position_row2)
                    break

                over = vcf.loc[position_row2]['POS']+len(vcf.loc[position_row2]['REF'])-1 + overlap_extra
                #print('minus=',window_sum,over,window_sum - over)
                if over > window_sum: #or 
                    window_sum = over 

                if window_sum - vcf.loc[position_row]['POS']>1000:
                    print('Warning interval contain large indel','interval ignored',[position_row,position_row2],'Indel lenght =',window_sum-vcf.loc[position_row]['POS'])#for filter large variant
                    if log_file:
                        log_file.write('Warning interval contain large indel, interval ignored ' + str([position_row,position_row2])+' Indel lenght = ' + str(window_sum-vcf.loc[position_row]['POS']) +'\n')
                    break

                #vcf_slice = vcf[(vcf.loc[position_row]['POS'] <= vcf['POS']) & (vcf['POS'] <= vcf.loc[position_row2]['POS']) & (vcf['#CHROM']==contig)]#!window_sum
                if type_merge == 'soft':
                    vcf_slice = vcf[(vcf.loc[position_row]['POS'] <= vcf['POS']) & (vcf['POS'] <= vcf.loc[position_row2]['POS']) & (vcf['#CHROM']==contig)]#!window_sum
                    #print(vcf_slice)
                    if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                        print('Warning interval contain indeterminate variant,interval ignored',[position_row,position_row2])
                        if log_file:
                            log_file.write('Warning interval contain indeterminate variant,interval ignored ' + str([position_row,position_row2])+'\n')
                        break

                interval = [vcf.loc[position_row]['POS'],window_sum]
                #interval2 = [vcf.loc[position_row]['POS'],vcf.loc[position_row2]['POS']]

                interval2 = [contig,vcf.loc[position_row]['POS'],vcf.loc[position_row2]['POS']]
                interval_dict[contig] = [vcf.loc[position_row]['POS'],vcf.loc[position_row2]['POS']]
            else:
                break
        
        if interval:
            interval_noexact.append(interval)
            
        if interval2:
            #interval_exact.append(interval2)
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
    pos_end = vcf['POS'] + vcf['REF'].apply(len) - 1

    for index,bedline in intervals_alignment.iterrows():
        if bedline['name'] in vcf.columns:
            if (bedline['start_position_ref']!=0 or bedline['end_position_ref']!=0) and (bedline['start_position_alt']!=0 or bedline['end_position_alt']!=0):
                vcf.loc[(vcf['POS']>=bedline['start_position_ref']) & (vcf['POS']<=bedline['end_position_ref']) & (vcf[bedline['name']]=='.'),bedline['name']] = 0#maybe bug with right unknown pos
                sum_0+=1
            else:
                sum_NA+=1
        else:
            print('Error sample not in vcf')
            
    return vcf



if __name__ == "__main__":


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

    #pd.set_option('display.max_columns', 10)


    local = False
    test = True
    local_windows = False
    info_just_indel = False


    if not(local or test or local_windows):
        intervals_path = parser.parse_args().interval.split(',')
        print('intervals_path=',intervals_path)
        ref_assemble_name = os.path.basename(file_fasta)[:-4]

    if local:
        directory = os.getcwd()

        file_vcf = '/home/strain5/Desktop/content/exp_super2_4/vcf_out/merged_exp_super2_4.vcf'
        
        file_fasta = '/home/strain5/Desktop/content/piplines/genomics_pipline_supply/AmesAncestor_GCF_0000084451.fna'
        file_gbk = '/home/strain5/Desktop/content/piplines/genomics_pipline_supply/' + 'AmesAncestor_GCF_000008445.1.gbk'
        work_dir = '/home/strain5/Desktop/content/exp_super2_4/vcf_out/'

        out_dir = '/home/strain5/Desktop/content/piplines/genomics_pipline_test/exp1/'
        out_file = out_dir + 'vcf_merged_777.vcf'
        log_file_path = out_dir + '/' + 'log_test2.txt'


        intervals_path = [work_dir+'/'+i for i in os.listdir(work_dir) if '.bed' in i]
        #intervals_path = ['/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_0.bed', '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_1.bed']

        ref_assemble_name = os.path.basename(file_fasta)[:-4]
    
    if local_windows:
        directory = os.getcwd()

        file_vcf = 'C:\\Users\\Grin\\Desktop\\remote_work\\exp_super2_4\\vcf_out\\merged_exp_super2_4.vcf'
        
        file_fasta = 'C:\\Users\\Grin\\Desktop\\remote_work\\piplines\\genomics_pipline_supply\\AmesAncestor_GCF_0000084451.fna'
        file_gbk = 'C:\\Users\\Grin\\Desktop\\remote_work\\piplines\\genomics_pipline_supply\\' + 'AmesAncestor_GCF_000008445.1.gbk'
        work_dir = 'C:\\Users\\Grin\\Desktop\\remote_work\\exp_super2_4\\vcf_out\\'

        out_dir = 'C:\\Users\\Grin\\Desktop\\remote_work\\piplines\\genomics_pipline_test\\exp2\\'
        out_file = out_dir + 'vcf_merged_777.vcf'
        log_file_path = out_dir + 'log_test111.txt'


        intervals_path = [work_dir+i for i in os.listdir(work_dir) if '.bed' in i]
        #intervals_path = ['/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_0.bed', '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA2_test/vcf_out/exp_A2_group_1.bed']

        ref_assemble_name = os.path.basename(file_fasta)[:-4]

    if test:
        directory = os.getcwd()

        file_fasta = os.path.join(directory,'test','test_merger.fna')
        file_vcf = os.path.join(directory,'test','test_merger_alignment_checker.vcf')
        work_dir = os.path.join(directory,'test')
        log_file_path = os.path.join(directory,'test','log.txt')

        out_dir = work_dir
        out_file = os.path.join(directory,'test','test_merger.vcf')
        intervals_path = [os.path.join(work_dir,i) for i in os.listdir(work_dir) if i.endswith('.bed')]

        file_gbk = os.path.join(os.path.split(directory)[0],'genomics_pipline_supply','AmesAncestor_GCF_000008445.1.gbk')

        #ref_assemble_name = 'A'


    log_file = open(log_file_path,'a')


    fasta = pipeline_base.SequenceFasta(file_fasta)
    fasta.seq_process(strip=True)
    sequence = ''.join(fasta.seq_lst)

    time_block_start = time.time()
    time_initial = time.time()


    header,vcf_head,type_head = pipeline_base.vcf_head_process(file_vcf)

    vcf = pd.read_csv(file_vcf,sep='\t',header = header,dtype=type_head)
    #vcf.index =vcf['POS'].values
    #vcf.index = vcf.index.astype(str) + '_' + vcf['#CHROM'].values.astype(str) + '_' +  vcf['POS'].values.astype(str)
    vcf.index = vcf['#CHROM'].astype(str) + '_' + vcf['POS'].astype(str) + '_' + vcf.index.astype(str)

    log_file = open(log_file_path,'a')
    log_file.write('\n'+'Start vcf_corrector'+'\n\n')
    print('\n'+'Start intervals_concat'+'\n\n')
    #intervals_alignment = intervals_concat(intervals_path)
    intervals_alignment = intervals_concat_bed(intervals_path)

    print('\n'+'Start vcf_corrector'+'\n\n')
    vcf_correct_bed = vcf_corrector_bed(vcf.copy(),intervals_alignment)
    timecheck('vcf_corrector_bed_after')
    

    log_file = open(log_file_path,'a')
    log_file.write('\n'+'Start definer_overlap_window'+'\n\n')
    print('\n'+'Start definer_overlap_window'+'\n\n')
    interval_exact,interval_noexact = definer_overlap_window(vcf_correct_bed,overlap_extra=0,log_file=log_file)
    #print('interval_noexact =',interval_noexact,'interval_exact=',interval_exact)

    vcf_correct_bed.to_csv(os.path.join(out_dir,'vcf_correct_bed.vcf'),sep='\t',index=False)

    log_file = open(log_file_path,'a')
    log_file.write('\n'+'Start merge_window'+'\n\n')
    print('\n'+'Start merge_window'+'\n\n')
    vcf_merged = merge_window(interval_exact,vcf_correct_bed.copy(),log_file=log_file,ref_sequence=sequence,fullcheck=False,intervals_alignment_bool=True,intervals_alignment=intervals_alignment)

    #vcf_merged = merge_window([['A',241, 251]],vcf_correct_bed.copy(),log_file=log_file,ref_sequence=sequence,fullcheck=False,intervals_alignment_bool=True,intervals_alignment=intervals_alignment)

    find_locus, find_source ,find_source_real = pipeline_base.contig_finder_gbk(file_gbk)
    vcf_merged = pipeline_base.position_editer(vcf_merged.copy(),find_locus,find_source)

    vcf_merged.to_csv(out_file,sep='\t',index=False)

    log_file.close()
    timecheck('end')
    print('end')