import os

import pandas as pd

print('start')

pd.set_option('display.max_columns', 20)

file_intervals1 = '/home/strain4/Desktop/genomics_pipline/interval1.txt'
file_intervals2 = '/home/strain4/Desktop/genomics_pipline/interval2.txt'

intervals_path = [i for i in os.listdir() if 'interval' in i]


file_fasta = '/home/strain4/Desktop/genomics_pipline/test_mini_vcf.fna'
file_fasta = '/home/strain4/Desktop/genomics_pipline/test_merger.fna'
file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger.vcf'
file_vcf = '/home/strain4/Desktop/genomics_pipline/test_merger_alignment_checker2.vcf'


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

    for pos_start,pos_end in intervals:
        print(pos_start,pos_end)
    
        vcf_slice = vcf[(pos_start <= vcf['POS']) & (vcf['POS'] <= pos_end) & (vcf['#CHROM']==contig)]

        vcf_slice_end = vcf[(vcf['POS'] > pos_end) & (vcf['#CHROM']==contig)]
        position_outside = int(vcf_slice_end.iloc[0]['POS'])

        sample_dict = {i:'' for i in vcf_slice.iloc[:,9:]}
        sample_bin_dict = {i:None for i in vcf_slice.iloc[:,9:]}

        position_first = int(vcf_slice.iloc[0]['POS'])
        position_start = int(vcf_slice.iloc[0]['POS'])
        position_last = int(vcf_slice.iloc[-1]['POS'])

        max_over = max([len(ref) + pos -1 for ref,pos in zip(vcf_slice['REF'],vcf_slice.index)])

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
        for position_row in vcf_slice.index:
            variant = [vcf_slice.loc[position_row]['REF']] + vcf_slice.loc[position_row]['ALT'].split(',')
            print(variant)
            position_inter = int(vcf_slice.loc[position_row]['POS'])
            #vcf.loc[position_row,'REF'] = ref_mod #maybe need
            for sample in sample_dict:

                if vcf_slice.loc[position_row][sample] == '.':
                    sample_dict[sample]
                    sample_will_deleted.add(sample)
                    #print(sample_dict,vcf_slice.loc[position_row][sample])
                else:
                    #print(variant)
                    sample_dict[sample] += sequence[position_start:position_inter-1] + variant[int(vcf_slice.loc[position_row][sample])]
            position_start = position_inter

        #for redefine bin in sample
        bin_value = 1
        used = []
        flag = False
        sample_lst = []

        for sample in sample_will_deleted:
            del sample_dict[sample]
            sample_bin_dict[sample] = '.'


        for sample in sample_dict:
            if all(vcf_slice[sample].astype(int) == 0): #!need check . values        
                sample_bin_dict[sample] = 0
            else:
                if sample not in used:
                    for sample_compared in sample_dict:
                        if vcf_slice.iloc[:,9:][sample_compared].equals(vcf_slice.iloc[:,9:][sample]) and sample_compared not in used:
                            sample_bin_dict[sample_compared] = bin_value
                            used += [sample_compared]
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

        for_merged = pd.Series(vcf_variant_dict)
        vcf.drop(vcf_slice.index,inplace=True)
        #vcf_merged = vcf.append(for_merged,ignore_index=True)
        vcf = vcf.append(for_merged,ignore_index=True)

    vcf.sort_values(by=['#CHROM','POS'],ascending=[False,True],inplace=True)

    return vcf

#vcf_merged = merge_window([[100,109]],contig,vcf.copy())
#vcf_merged = merge_window(294715,294722,'NC_007530')

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
            #print('loop1=','position_row=',position_row,'window_sum=',window_sum,'position_row2=',position_row2,'over=',over)            

            if position_row2 <= window_sum:

                #if any(vcf.loc[position_row2][9:] == '.'):
                #    continue

                over = position_row2+len(vcf.loc[position_row2]['REF'])-1
                #print('loop2=','position_row=',position_row,'window_sum=',window_sum,'position_row2=',position_row2,'over=',over)            
                if over > window_sum:
                    window_sum = over

                #interval = [position_row,window_sum]
                #interval2 = [position_row,position_row2]

                vcf_slice = vcf[(position_row <= vcf['POS']) & (vcf['POS'] <= position_row2) & (vcf['#CHROM']==contig)]
                if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                    #print('EEEE3',vcf_slice)
                    break
                interval = [position_row,window_sum]
                interval2 = [position_row,position_row2]                

            else:
                #vcf_slice = vcf[(position_row <= vcf['POS']) & (vcf['POS'] <= window_sum) & (vcf['#CHROM']==contig)]
                #if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
                #    print('position_row=',position_row,'over=',over,'window_sum=',window_sum)
                #    print('EEEE2',vcf_slice)
                #    interval = []
                #    interval2 = []
                break

        if over - position_row>1000:#for filter large variant
            break

        #if interval or interval2:
        #    vcf_slice = vcf[(interval2[0] <= vcf['POS']) & (vcf['POS'] <= interval2[1]) & (vcf['#CHROM']==contig)]
        #    if any((vcf_slice.iloc[:,9:].isin(['.']).any()) & (~(vcf_slice.iloc[:,9:].isin(['.']).all()))):
        #        print('EEEE4',vcf_slice)
        #        break

        #print('interval=',interval)
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

'''for interval in interval_noexact:
    intervals_alignment_slice = intervals_alignment[(intervals_alignment['A_start'] < interval[0]) & (interval[1] < intervals_alignment['A_end'])]
    
    #if not all(vcf_slice[sample].isin(['.'])) and not all(~vcf_slice[sample].isin(['.']))
    #(~(intervals_alignment_slice == 0)).all(axis=1)
    #(intervals_alignment_slice == 0).all(axis=1)
    (intervals_alignment_slice == 0 | intervals_alignment_slice.isna())
    (intervals_alignment_slice == 0 | intervals_alignment_slice.isna()).all(axis=1)
    if any((intervals_alignment_slice == 0).any(axis=1)):
    #if not all((~(intervals_alignment_slice == 0)).all(axis=1)) and not all((intervals_alignment_slice == 0).all(axis=1)):
        print('cool')
    else:
        print('no')'''



    #print(intervals_alignment_slice)

#print(interval_exact)
print('interval_noexact =',interval_noexact)


#vcf_merged = merge_window(interval_noexact,vcf.copy(),contig)
vcf_merged = merge_window([[100,103]],vcf.copy(),contig)

#vcf_merged = merge_window([[81,85]],vcf.copy(),contig)

#vcf_merged = merge_window([[100,109]],vcf.copy(),contig,)

#vcf_merged.to_csv('/home/strain4/Desktop/test_split/test_split_m2.vcf',sep='\t')


print('end')