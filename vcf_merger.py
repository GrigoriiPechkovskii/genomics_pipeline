import os

import pandas as pd

print('start')

pd.set_option('display.max_columns', 20)

file_intervals1 = '/home/strain4/Desktop/genomics_pipline/interval1.txt'
file_intervals2 = '/home/strain4/Desktop/genomics_pipline/interval3.txt'

intervals_path = [i for i in os.listdir() if 'interval' in i]


file_fasta = '/home/strain4/Desktop/genomics_pipline/test_mini_vcf.fna'
file_fasta = '/home/strain4/Desktop/genomics_pipline/test_merger.fna'
file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger.vcf'
file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker4.vcf'
file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker3.vcf'


contig = 'A'
header = 5
header = 0
#file_fasta = '/home/strain4/Desktop/test_split/GCF_000008445.1_ASM844v1_genomic.fna'
#file_vcf = '/home/strain4/Desktop/test_split/merged.vcf'
#contig = 'NC_007530'
#header = 15

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
vcf = pd.read_csv(file_vcf,sep='\t',header = header)#!15
vcf.index =vcf['POS'].values

sequence = fasta.seq_lst[0]

#def merge_window(pos_start,pos_end,contig,vcf,ref_sequence=sequence):
def merge_window(intervals,vcf,contig,ref_sequence=sequence):
    '''Функция которая обьединяет заданое окно в vcf файле по позициям'''
    global vcf_slice
    intervals_used = []
    for pos_start,pos_end in intervals:

        

        if any([pos_start_used <= pos_start <= pos_end_used for pos_start_used,pos_end_used in intervals_used]):
            print('Warning this interval used!',[pos_start,pos_end])
            continue

        intervals_used.append([pos_start,pos_end])
        print('INTERVALS===',[pos_start,pos_end])


        print(pos_start,pos_end)
    
        vcf_slice = vcf[(pos_start <= vcf['POS']) & (vcf['POS'] <= pos_end) & (vcf['#CHROM']==contig)]

        vcf_slice_end = vcf[(vcf['POS'] > pos_end) & (vcf['#CHROM']==contig)]
        position_outside = int(vcf_slice_end.iloc[0]['POS'])

        #sample_dict = {i:'' for i in vcf_slice.iloc[:,9:]}

        sample_dict_mem = {i:0 for i in vcf_slice.iloc[:,9:]}

        sample_bin_dict = {i:None for i in vcf_slice.iloc[:,9:]}

        position_first = int(vcf_slice.iloc[0]['POS'])
        position_start = int(vcf_slice.iloc[0]['POS'])
        position_last = int(vcf_slice.iloc[-1]['POS'])


        max_over = max([len(ref) + pos - 1 for ref,pos in zip(vcf_slice['REF'],vcf_slice.index)])


        sample_dict = {i:list(sequence[position_first-1:max_over]) for i in vcf_slice.iloc[:,9:]}
        sample_dict_mem = {i:position_first for i in vcf_slice.iloc[:,9:]}

        print('position_first=',position_first,'max_over=',max_over,'sample_dict=',sample_dict,vcf_slice.index)

        position_last = max_over

        ref_mod = sequence[position_first-1:position_last]

        for position_row in vcf_slice.index[:-1]:
            #print('w =',position_row + len(vcf_slice.loc[position_row]['REF']),position_outside)
            if position_row + len(vcf_slice.loc[position_row]['REF']) - position_outside > 0:
                print('Warning overlay window maybe wrong result,','position =',position_row) 
                #print(position_row + len(vcf_slice.loc[position_row]['REF']) - position_last)        

        flag = False
        if  any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
            print('EEEE',vcf_slice)

        for sample in sample_dict:
            if  not all(vcf_slice[sample].isin(['.'])) and not all(~vcf_slice[sample].isin(['.'])):#!
                #print('Warning indeterminate genotype of the sample')
                #assert False, "'Warning indeterminate genotype of the sample'"
                flag = True
        if flag:
            continue                

        variant = []
        sample_will_deleted = set()
        print('st=',sample_dict)
        for position_row in vcf_slice.index:
            variant = [vcf_slice.loc[position_row]['REF']] + vcf_slice.loc[position_row]['ALT'].split(',')
            position_inter = int(vcf_slice.loc[position_row]['POS'])
            #vcf.loc[position_row,'REF'] = ref_mod #maybe need
            for sample in sample_dict:

                if vcf_slice.loc[position_row][sample] == '.':
                    sample_dict[sample]
                    sample_will_deleted.add(sample)
                else:
                    #if sample_dict_mem[sample] > position_row:
                        #print(position_inter)
                    #if sample_dict_mem[sample] < max_over and int(vcf_slice.loc[position_row][sample]) != 0 or True:
                    if int(vcf_slice.loc[position_row][sample]) != 0:
                        #if vcf_slice.loc[position_row][sample] == '0' or vcf_slice.loc[position_row][sample] == 0:# !type
                            #continue
                        #else:
                            #n = max_over - position_row
                        n = position_row - position_first
                        m = len(variant[int(vcf_slice.loc[position_row][sample])])
                        r_l = len(vcf_slice.loc[position_row]['REF'])
                        
                        print('m=', len(variant[int(vcf_slice.loc[position_row][sample])]),'r_l=',r_l,'n=',n)
                        if len(variant[int(vcf_slice.loc[position_row][sample])]) == r_l:#!

                            sample_dict[sample][n:n+r_l] =  [list(variant[int(vcf_slice.loc[position_row][sample])])]

                        elif len(variant[int(vcf_slice.loc[position_row][sample])]) != r_l:
                            print('+++',list(variant[int(vcf_slice.loc[position_row][sample])]) + [''] * (r_l - m ))
                            #sample_dict[sample][n:n+r_l] = list(variant[int(vcf_slice.loc[position_row][sample])]) + [''] * (r_l - m)
                            
                            sample_dict[sample][n:n+r_l] = [list(variant[int(vcf_slice.loc[position_row][sample])])] + [''] * (r_l - 1)
                            

                        sample_dict_mem[sample] += len(vcf_slice.loc[position_row]['REF'])
                            
                        print('f=',sample,position_row,sample_dict,'\n',sample_dict_mem,variant[int(vcf_slice.loc[position_row][sample])],vcf_slice.loc[position_row][sample],'\n',position_start,position_inter)

                    else:
                        print('sample_dict_mem[sample] >= max_over',sample_dict_mem[sample],max_over,'alt=',int(vcf_slice.loc[position_row][sample]),'sample=',sample)
                        continue
            print('end circle',sample_dict,'\n',sample_dict_mem,variant[int(vcf_slice.loc[position_row][sample])],vcf_slice.loc[position_row][sample],'\n',position_start,position_inter)

            position_start = position_inter
        
        #for redefine bin in sample
        bin_value = 1
        used = []
        used_bin = []
        flag = False
        sample_lst = []

        for sample in sample_will_deleted:
            del sample_dict[sample]
            sample_bin_dict[sample] = '.'

        print('3',sample_dict)


        for sample in sample_dict:
            sample_dict[sample] = [i if type(i) == str else ''.join(i) for i in sample_dict[sample]]
        sample_dict = {i:''.join(j) for i,j in sample_dict.items()}

        print('4',sample_dict)

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
        vcf_variant_dict = {'#CHROM':contig,'POS':position_first,'ID':'.','REF':ref_mod,
                    'ALT':alt,'QUAL':40,'FILTER':'PASS','INFO':info,'FORMAT':'GT'}
        vcf_variant_dict.update(sample_bin_dict)
        #print('vcf_variant_dict=',vcf_variant_dict)

        for_merged = pd.Series(vcf_variant_dict,name=position_first)
        vcf.drop(vcf_slice.index,inplace=True)
        #vcf_merged = vcf.append(for_merged,ignore_index=True)
        vcf = vcf.append(for_merged,ignore_index=False)#!True

    vcf.sort_values(by=['#CHROM','POS'],ascending=[False,True],inplace=True)

    return vcf


def definer_overlap_window(vcf):
    over = 0
    interval_noexact = []
    interval_exact = []
    pos = vcf.iloc[0]["POS"]
    for num,position_row in enumerate(vcf.loc[pos:].index):
        interval = []
        interval2 = []
        window_sum = position_row + len(vcf.loc[position_row]['REF'])-1 
        for position_row2 in vcf.index[num+1:]:
            if position_row2 <= window_sum:
                over = position_row2+len(vcf.loc[position_row2]['REF'])-1
                if over > window_sum:
                    window_sum = over

                vcf_slice = vcf[(position_row <= vcf['POS']) & (vcf['POS'] <= position_row2) & (vcf['#CHROM']==contig)]
                if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                    break
                interval = [position_row,window_sum]
                interval2 = [position_row,position_row2]
            else:
                break

        if over - position_row>1000:#for filter large variant
            break
        
        if interval:
            interval_noexact.append(interval)
            
        if interval2:
            interval_exact.append(interval2)

    return interval_exact,interval_noexact #interval_exact not work

def intervals_concat(intervals_path):
    '''intervals_path = list of path file interval'''
    intervals = []
    for interval in intervals_path:
        intervals.append(pd.read_csv(interval,sep='\t'))
    intervals_full = pd.concat(intervals,sort=False)

    return intervals_full

intervals_alignment = intervals_concat(intervals_path)

interval_exact,interval_noexact = definer_overlap_window(vcf)

print('interval_noexact =',interval_noexact,'interval_exact=',interval_exact)


#vcf_merged = merge_window(interval_exact,vcf.copy(),contig)
#vcf_merged = merge_window([[100,107]],vcf.copy(),contig)
vcf_merged = merge_window([[161,170]],vcf.copy(),contig)
#vcf_merged.to_csv('/home/strain4/Desktop/test_split/test_split_m2.vcf',sep='\t')


print('end')