#!/usr/bin/python3

#By Grigorii Pechkovskii
'''
Check ext with blank alt - fixed, its maybe inversion, out with warning string
Loging
Reference genome must be flexible
N letters
write test position between original fasta and alignment

#!bug in l Warning Exception snr start with snp 994529 252142 252136 {'T', 'G'}
#!bug reapet can over position each other
#mod test please
'''
print('start')

import re
import os
import sys
import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-x', '--xmfa',action='store', help='File xmfa')
parser.add_argument('-r', '--ref',action='store', help='Reference')
parser.add_argument('-d', '--dir',action='store', help='Work directory')
parser.add_argument('-o', '--out-dir',action='store',default=os.getcwd(), help='Out directory')
parser.add_argument('-i', '--index-type',action='store',default='mauve', help="Index type 'mauve' or 'parsnp'")
parser.add_argument('-g', '--gbk-file',action='store',default='mauve', help='File gbk')
parser.add_argument('-n', '--name-vcf',action='store',default='test_mini.vcf', help='File gbk')

REF = parser.parse_args().ref
directory_file_xmfa = parser.parse_args().xmfa
directory  = parser.parse_args().dir
directory_out  = parser.parse_args().out_dir
file_gbk = parser.parse_args().gbk_file
index_type = parser.parse_args().index_type
name_vcf = parser.parse_args().name_vcf

if not os.access(directory_out,os.F_OK):
    os.mkdir(directory_out)

if True:
    directory = os.getcwd()
    directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/test_mini.xmfa'
    #directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/parsnp.xmfa'
    #directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/mauve_out1.xmfa'
    #directory_file_xmfa = '/home/strain4/Desktop/content/GI_AAJ/GI_AAJ_out1/out1'
    #directory_file_xmfa = '/home/strain4/Desktop/content/GI_AAF/choise3_out2/parsnp.xmfa'
    #directory_file_xmfa = '/home/strain4/Desktop/piplene_mauve/exp1/exp_0/group0'
    directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/exp_1_group_1'


    REF = 'AmesAncestor_GCF_000008445.1'#test_mini
    #REF ='AmesAncestor_GCF_0000084451'
    REF = 'GCF_000008445.1_ASM844v1_genomic'#GI_AAJ_out1
    #REF = 'Ames_Ancestor_ref_GCF_000008445.1_ASM844v1_genomic.fna'#parsnp.xmfa
    directory_out = os.getcwd()


    index_type = 'mauve'
    #name_vcf_simple = 'test_sim.vcf'

    file_gbk = directory + '/AmesAncestor_GCF_000008445.1.gbk'
    name_vcf = 'test3.vcf'

sort = True
delete_ref = True


if index_type == 'parsnp':    
    pos_vcf = 1
    pos_minus = 1
elif index_type == 'mauve':
    pos_vcf = 0
    pos_minus = 2

def get_index(directory_file_xmfa,index_type):
    '''Iter on xmfa header(mauve format) with #
        and return dict with id(key) and name genome(values)
    '''
    file_xmfa = open(directory_file_xmfa)
    id_nameseq_dict = dict()

    if index_type == 'mauve':
        for file_line in file_xmfa:
            if '#' in file_line:
                find = re.search(r'Sequence(\d*)File\t(.*)\.',file_line)
            else:
                break
            if find != None:
                id_nameseq_dict[find.group(1)] = os.path.basename(find.group(2))

    if index_type == 'parsnp':
        find_id_lst = []
        find_genome_lst = []
        for file_line in file_xmfa:
            if '#' in file_line:
                #print(id_nameseq_dict)
                find_id = re.search(r'SequenceIndex\s(\d*)',file_line)
                find_genome = re.search(r'SequenceFile\s(.*)\.',file_line)#!\w replace .
            else:                
                id_nameseq_dict = dict(zip(find_id_lst,find_genome_lst))
                break

            if find_id != None:
                find_id_lst.append(find_id.group(1))
            if find_genome != None:
                find_genome_lst.append(find_genome.group(1))

            #if '>' in file_line:

                    
    file_xmfa.close()
    return id_nameseq_dict


id_nameseq_dict = get_index(directory_file_xmfa,index_type)
id_nameseq_dict_val = list(id_nameseq_dict.values())

def head_vcf(name_seq:list):
    '''Function for make vcf header,
       get name_seq: list of name from xmfa header with strict order
    '''
    columns_name =['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    columns_name += name_seq
    tag = '##fileformat=VCFv4.2' + '\n' 
    tag += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n'
    tag += '##FILTER=<ID=PASS,Description="All filters passed">' + '\n'
    tag += '##INFO=<ID=SNP,Number=0,Type=Flag,Description="Indicates that the variant is an SNP.">' + '\n'
    tag += '##INFO=<ID=SNR,Number=0,Type=Flag,Description="Indicates that the variant is an SNR.">' + '\n'
    tag += '##INFO=<ID=REPEAT,Number=0,Type=Flag,Description="Indicates that the variant is an REPEAT.">' + '\n' 
    tag += '##INFO=<ID=ComplexIndel,Number=0,Type=Flag,Description="Indicates that the variant is an ComplexIndel.">' + '\n' 
    tag += '##INFO=<ID=NANinfo,Number=0,Type=Flag,Description="Indicates that the variant is an NANinfo.">' + '\n' 
    tag += '##INFO=<ID=substitution,Number=0,Type=Flag,Description="Indicates that the variant is an substitution.">' + '\n' 
    tag += '##INFO=<ID=maybe_inversion,Number=0,Type=Flag,Description="Indicates that the variant is an maybe_inversion.">' + '\n'   
    head = '\t'.join(columns_name) + '\n'
    with open(directory_out +'/'+ name_vcf,'w') as file_vcf:
        file_vcf.write(tag)
        file_vcf.write(head)
    return tag
    #print(head)

tag = head_vcf(id_nameseq_dict_val)

def single_aln_generator(directory_file_xmfa):
    '''Generator for xmfa,
       get a directory_file_xmfa
       yield separate with = aln '''
    title_seq = []
    seq_seq = []
    single_aln = ''
    notfirst = False
    tmp = ''
    file_xmfa = open(directory_file_xmfa)    
    for line in file_xmfa:
        if '#' not in line and '=' not in line: #dont want make check for all string
            single_aln += line                
            if '>' in line:
                title_seq.append(line.strip())
                if notfirst:
                    seq_seq.append(tmp)
                    tmp = ''                
            else:
                notfirst = True
                tmp += line.strip().upper() #!up register 

        if '=' in line:
            seq_seq.append(tmp)
            tmp = ''
            notfirst = False
            yield title_seq , seq_seq
            title_seq = []
            seq_seq = []

    file_xmfa.close()


def parser_title(title_seq:list):
    '''Parsing > line
       get title_seq(> line)
       return positioin start - end sequence
       and name_idseq - name sequece from header name and id in > line
    '''
    position_every = []
    position_every_dict = dict()
    position_every_num_dict = dict()
    name_every = []
    id_strand = dict()
    name_strand = dict()
    for title in title_seq:

        id_search = re.search(r'(\d*):',title)
        id_seq = id_search.group(1)

        #!!!MAUVE
        #name_search = re.search(r'(/.*)\.',title)
        #name = name_search.group(1)
        #name = os.path.basename(name)
        #name_every.append(name)

        position_search = re.search(r'\d*:(\d*)-(\d*)\b',title)
        position = [position_search.group(1),position_search.group(2)]
        position_every += [position]

        #position_every_dict[name] = position

        position_every_num_dict[id_seq] = position

        strand_search = re.search(r'\s([\+-])\s',title_seq[0]).group(1)
        id_strand[id_seq] = strand_search


    #print(id_strand)
    name_seq_title = dict()
    #for sorting maybe!!
    for key,val in position_every_num_dict.items() :
        #print(id_nameseq_dict[key])
        name_seq_title[id_nameseq_dict[key]] = val#!id_nameseq_dict frome up namespace
        

    for key,val in id_strand.items():
        name_strand[id_nameseq_dict[key]] = val

    pos = list(name_seq_title.values())
    name_idseq = list(name_seq_title)
    strand = list(name_strand.values())
    #print(name_strand)
    #return name_every,position_every,position_every_dict,position_every_num_dict,pos,name_idseq
    return pos,name_idseq,strand


def join_dif(massive_str):
    massive_list = []
    for one_str in massive_str:
        massive_list.append(''.join(list(one_str)))
    return massive_list


def snr_diffinder2(sym_num,start_pos,sym_seq_lst,info):
    
    ref_pos = int(pos_set[0][0][0])

    info_flag = False

    join_dif(sym_seq_lst[1])

    sym_seq_lst[1:]
    alt_set_tmp = set()
    any(list((s==['-']*len(s) for s in sym_seq_lst[1:])))
    
    list((alt_set_tmp.update(set(s)) for s in sym_seq_lst[1:]))
    ref_seq = seq_seq[0]
    #if len(sym_seq_lst) == 2 and '-' in sym_seq_lst[1] and all([len(set(s))==2 for s in sym_seq_lst[1:]]): #len(set(sym_seq_lst[1])) == 2:# any('-' in s for s in sym_seq_lst[1]):
    
    if any('-' in s for s in sym_seq_lst) and len(alt_set_tmp) == 2:
        #print(sym_seq_lst)
        if any([s == tuple('-')*len(sym_seq_lst[1:]) for s in zip(*sym_seq_lst[1:])]):
            print(sym_seq_lst)

        info  = 'SNR'
        uniq_let_lst = [i for i in sym_seq_lst[1] if i!='-']
        uniq_let_set = set(uniq_let_lst)
        uniq_let = uniq_let_lst[0]


        uniq_lst_right=[]
        for right_ind in range(sym_num,len(ref_seq)-1):
            if all([(seq[right_ind]==uniq_let or seq[right_ind]=='-' )for seq in seq_seq]):#!!!ex if break in snp
                uniq_lst_right.append([uniq_let]*len(seq_seq))
            else:   
                break
        if len(set([seq[right_ind] for seq in seq_seq]))>1:
                    print('r Warning Exception snr ended with snp',start_pos)


        uniq_lst_left=[]
        for left_ind in range(sym_num-2,0,-1):
            if all([(seq[left_ind]==uniq_let or seq[left_ind]=='-' )for seq in seq_seq]):#!!!ex if break in snp
                uniq_lst_left.append([uniq_let]*len(seq_seq))
            else:
                break
        left_end = [ref_seq[left_ind]]*len(seq_seq)
        left_end_pos = left_ind
        if len(set([seq[left_ind] for seq in seq_seq]))>1:
                    print('l Warning Exception snr start with snp',start_pos,sym_num-2,left_ind,set([seq[left_ind] for seq in seq_seq]))


        sym_seq_lst = [left_end] +  uniq_lst_left + sym_seq_lst[1:] + uniq_lst_right
        start_pos = left_end_pos + ref_pos


    return start_pos,sym_seq_lst,info

def tandem_check(list_nuc_strings:list):
    '''Get list of strings
       return list of min tandam repeat'''
    list_tandem = []

    for nuc_string in set(list_nuc_strings):

        result = re.findall(r'(\w+?)\1+', nuc_string)
        if '-' not in nuc_string:            

            if len(result) != 0:
                for s in result:
                    if nuc_string.replace(s,'') == '':
                        list_tandem.append(s)
                    else:
                        list_tandem.append(nuc_string)

            

            elif len(result) == 0:
                list_tandem.append(nuc_string)

    if len(set(list_tandem)) == 1:

        #print('tandem',list_tandem[0])
        return list_tandem[0]
    else:
        #print('Warning in func tandem_check not 1 alt tandem repeat',list_tandem)
        pass
        #return list_tandem[0]#!return complex uniq_let
            

def repeat_diffinder(sym_num,start_pos,sym_seq_lst,info,ref_pos,ref_seq):
    
    #ref_pos = int(pos_set[0][0][0])

    info_flag = False
    pos_shift = -1
    #join_dif(sym_seq_lst[1])

    alt_set_tmp = set()
    any(list((s==['-']*len(s) for s in sym_seq_lst[1:])))
    
    list((alt_set_tmp.update(set(s)) for s in sym_seq_lst[1:]))
    #ref_seq = seq_seq[0]
    #print(start_pos,sym_seq_lst)
    #print(start_pos,[s  for s in sym_seq_lst[1:]])
    if any([s == tuple('-')*len(sym_seq_lst[1:]) for s in zip(*sym_seq_lst[1:])]):
        uniq_let = tandem_check(join_dif(zip(*sym_seq_lst[1:])))
        if uniq_let == None:
            info = 'ComplexIndel'
            return start_pos,sym_seq_lst,info
        
        uniq_lst_right=[]
        #print(sym_num,start_pos)
        #for right_ind_start,right_ind_end in zip(range(sym_num+1+len(uniq_let),len(ref_seq),len(uniq_let)),range(sym_num+len(uniq_let)*2,len(ref_seq)-1,len(uniq_let))):
        log_bug = open('log_bug','a')
        print('start_pos',start_pos,'sym_num=',sym_num,'uniq_let=',uniq_let,file=log_bug) 
        for right_ind_start,right_ind_end in zip(range(sym_num+len(sym_seq_lst[1:]),len(ref_seq)-1,len(uniq_let)),range(sym_num+len(sym_seq_lst[1:])+len(uniq_let),len(ref_seq)-1,len(uniq_let))):
            
            print(right_ind_start,right_ind_end,seq_seq[0][right_ind_start:right_ind_end],
                all([(seq[right_ind_start:right_ind_end]==uniq_let or seq[right_ind_start:right_ind_end]=='-'*len(uniq_let) )for seq in seq_seq]),file=log_bug)

            if all([(seq[right_ind_start:right_ind_end]==uniq_let or seq[right_ind_start:right_ind_end]=='-'*len(uniq_let) )for seq in seq_seq]):#!!!ex if break in snp
                uniq_lst_right.append([uniq_let]*len(seq_seq))
                #info  = 'SNR'
                if len(uniq_let) == 1:
                    info = 'SNR'
                else:
                     info = 'REPEAT'
            else:
                if len(set([seq[right_ind_start] for seq in seq_seq]))>1:
                    print('r Warning Exception snr ended with snp',start_pos)
                break
        
        uniq_lst_left=[]
        #left_end = [ref_seq[sym_num-1]]*len(seq_seq)
        for left_ind_start,left_ind_end in zip(range(sym_num,0,-len(uniq_let)),range(sym_num-len(uniq_let),0,-len(uniq_let))):
        #for left_ind_start,left_ind_end in zip(range(sym_num,0,-len(sym_seq_lst[1:])),range(sym_num-len(uniq_let),0,-len(uniq_let))):
            #[(print(seq[left_ind_end:left_ind_start],uniq_let,'-'*len(uniq_let),seq[left_ind_end:left_ind_start]==uniq_let,seq[left_ind_end:left_ind_start]=='-'*len(uniq_let),all([(seq[left_ind_end:left_ind_start]==uniq_let or seq[left_ind_end:left_ind_start]=='-'*len(uniq_let) )for seq in seq_seq]))) for seq in seq_seq]
            if all([(seq[left_ind_end:left_ind_start]==uniq_let or seq[left_ind_end:left_ind_start]=='-'*len(uniq_let) )for seq in seq_seq]):#!!!ex if break in snp
                uniq_lst_left.append([uniq_let]*len(seq_seq))
                #info  = 'SNR'
                if len(uniq_let) == 1:
                    info = 'SNR'
                else:
                     info = 'REPEAT'
                pos_shift += len(uniq_let)
            else:
                pos_shift += 1

                left_end = [ref_seq[left_ind_start-1]]*len(seq_seq)
                #left_end_pos = left_ind_start-1

                if len(set([seq[left_ind_start-1] for seq in seq_seq]))>1:
                    print('l Warning Exception snr start with snp',start_pos,sym_num,left_ind_start,set([seq[left_ind_start-1] for seq in seq_seq]))
                break
        else:
            left_end = [ref_seq[sym_num-1]]*len(seq_seq)
            pos_shift = 0
        #left_end = [ref_seq[left_ind_start-2]]*len(seq_seq)
        #left_end_pos = left_ind_start-2
        print('[left_end=]',[left_end],'uniq_lst_left=',uniq_lst_left,'sym_seq_lst[1:]',sym_seq_lst[1:],
            'uniq_lst_right=',uniq_lst_right,
            '+++=',[left_end] +  uniq_lst_left + sym_seq_lst[1:] + uniq_lst_right,file=log_bug)
        sym_seq_lst = [left_end] +  uniq_lst_left + sym_seq_lst[1:] + uniq_lst_right
        start_pos = start_pos - pos_shift


    #print(start_pos,sym_seq_lst)
    #log_bug.close()
    return start_pos,sym_seq_lst,info

def diffinder(seq_seq,ref_pos,ref_seq,pos_vcf=pos_vcf,pos_minus=pos_minus):
    #req add position parametr if need break loop through equal symbol
    #start_pos - ref position last equal symbol + pos vcf
    #rep_pos - sym_num position first non equal symbol
    sym_seq_start = 'YYY'
    info = 'NANinfo'
    #sym_seq_end = ''
    sym_num_start = 0
    sym_seq_lst = list()

    #ref_pos = int(pos_set[0][0][0])
    #ref_seq = seq_seq[0]#!reference sequence for compare  
    #pos_vcf = 0
    start_pos = ref_pos
    first_flag = True
    ref_flag = True
    rep_pos =0
    for sym_num,sym_seq in enumerate(zip(*seq_seq)):
        if ref_seq[sym_num] != '-':
            pos_vcf += 1
        if len(set(sym_seq)) > 1:#input in loop if alternative exist
            sym_seq_lst += [list(sym_seq)]
            
            if ref_flag:
                ref_flag = False
                #start_pos = ref_pos + pos_vcf - 1
                rep_pos = sym_num

        if (len(set(sym_seq)) == 1 and len(sym_seq_lst) != 0) or sym_num == len(ref_seq)-1:#here break
            ref_flag = True
            if sym_num == len(ref_seq)-1 and (len(set(sym_seq)) == 1):#why
                print('WARNING break in diffinder')
                break


            if first_flag:# and seq_seq[1][0:10]:
                first_flag = False
                #start_pos = ref_pos + pos_vcf - 2#end pos??
                #start_pos = ref_pos
                #if sym_num==1:
                if sym_num != len(ref_seq)-1:                    
                    #start_pos = ref_pos    
                    if len(sym_seq_lst) == sym_num:
                        start_pos = ref_pos + pos_vcf - 2
                        print('Warning origin seq start with indel 1',start_pos)
                    #if any([seq_seq[i][0:sym_num] == '-'*sym_num for i in range(len(seq_seq))]):#S.startswith(str)
                        if len(sym_seq_lst)==1 and all(['-' not in s for s in sym_seq_lst]):#ex for first snp
                            start_pos = ref_pos
                            #sym_seq_lst = [list(sym_seq_start)] + sym_seq_lst
                        else:
                            start_pos = ref_pos
                            sym_seq_lst = sym_seq_lst + [list(list(zip(*seq_seq))[sym_num])]
                    else:
                        sym_seq_lst = [list(sym_seq_start)] + sym_seq_lst
                        #start_pos = ref_pos + pos_vcf - 2

                elif sym_num == len(ref_seq)-1:
                    #sym_seq_lst = [list(sym_seq_start)] + sym_seq_lst
                    start_pos = ref_pos

            #if first_flag and sym_num !=0: !!
            #elif not first_flag and 
            else:               
                sym_seq_lst = [list(sym_seq_start)] + sym_seq_lst
            start_pos,sym_seq_lst,info = repeat_diffinder(rep_pos,start_pos,sym_seq_lst,info,ref_pos,ref_seq)
            yield start_pos,sym_seq_lst,info

            info = 'NANinfo'

            sym_seq_lst = []
            start_pos = ref_pos + pos_vcf - 1
            rep_pos = sym_num
            sym_seq_start = sym_seq            

        if not (len(set(sym_seq)) > 1) and not ((len(set(sym_seq)) == 1 and len(sym_seq_lst) != 0) or sym_num == len(ref_seq)-1):
            ref_flag = True
            sym_seq_start = sym_seq
            start_pos = ref_pos + pos_vcf - 1
            rep_pos = sym_num

        
def dif_process(seq_lst,position,info='NANinfo'):
    variance_bin_dict = dict()
    #info = 'NANinfo'
    name_str_dict = {x: "NAN" for x in id_nameseq_dict_val}
    variance_lst = join_dif(list((zip(*seq_lst))))
    variance_set = set((zip(*seq_lst)))

    alt_variance_set  = set()
    alt_variance_dict = dict()
    #print(variance_lst)
    if len(variance_lst)>0:
        reference_variance = variance_lst[0] #define reference

    alt_num = 1
    alt_variance_dict[reference_variance] = '0'

    alt_variance_set.add(reference_variance)

    for variance in variance_lst[1:]:
        if variance not in alt_variance_set:

            alt_variance_set.add(variance)
            alt_variance_dict[variance] = str(alt_num)
            alt_num += 1

    for variance_n in range(len(variance_lst)):
        name_str_dict[pos_set[1][variance_n]] = variance_lst[variance_n]

    for name in name_str_dict:
        if name_str_dict[name] != 'NAN':
            name_str_dict[name] = alt_variance_dict[name_str_dict[name]]

    #print(any(('N' in s for s in list(alt_variance_dict.keys()))))
    #print(list(alt_variance_dict.keys()))
    
    alt_variance_dict_n = {}
    #for SNP    
    #if all(len(s) == 2 for s in list(alt_variance_dict.keys())) and all('-' not in s for s in list(alt_variance_dict.keys())):  
    if len(seq_lst)>=2 and all('-' not in s for s in seq_lst):  
        for alt_key,alt_val in alt_variance_dict.items():
            alt_variance_dict_n[alt_key[1:]] = alt_val
        position += 1
        if len(seq_lst)==2:
            info = 'SNP'
        else:
            info = 'substitution'

        return name_str_dict,list(alt_variance_dict_n),position,info
    
    #!attention with dict maybe bag
    dash_logic = True

    if dash_logic:

        alt_variance_dict_r = {j:i.replace('-','') for i,j in alt_variance_dict.items()}#rearrangement dict  
        alt_variance_dict_n = {i.replace('-',''):j for i,j in alt_variance_dict.items()}
        #print('1',alt_variance_dict)
        #print('r',alt_variance_dict_r)
        #print('n',alt_variance_dict_n)

        #for alt_key,alt_val in alt_variance_dict_r.items():
            
            #alt_variance_dict_n[alt_key] = alt_val.replace('-','')
            #alt_variance_dict_r[alt_val.replace('-','')] = alt_key

        if len(set(alt_variance_dict_r.values())) != len(list(alt_variance_dict_r.values())):           

            info = 'maybe_inversion'
            #print('WARNING maybe wrong alignment')

            return name_str_dict,list(alt_variance_dict_r.values()),position,info
        return name_str_dict,list(alt_variance_dict_n),position,info

    return name_str_dict,list(alt_variance_dict),position,info#!maybe not


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

find_locus, find_source ,find_source_real = contig_finder_gbk(file_gbk)

def contig_definder(position,find_locus,find_source): 
    ''' ''' 
    for locus,source in zip(find_locus,find_source):
        if (source[0] <= position <= source[1]):
            position_real = position - source[0]+ 1#!
            return locus,position_real

'''file_genome = 'AmesAncestor_GCF_0000084451.fna'
file_genome_opened = open(file_genome)
file_genome_read = ''
for line in file_genome_opened:
    if '>' not in line:
        file_genome_read += line.strip()
file_genome_opened.close()'''


def seq_reverse(seq_seq):

    #seq_rev = ''
    seq_seq_rev = []
    for seq in seq_seq:
        seq_rev = ''
        for sym in seq:
            if sym=='C':
                seq_rev += 'G'
            elif sym=='G':
                seq_rev += 'C'
            elif sym=='T':
                seq_rev+='A'
            elif sym=='A':
                seq_rev += 'T'
            elif sym=='-':
                seq_rev += '-'
        seq_rev = seq_rev[::-1]
        seq_seq_rev += [seq_rev]
    return seq_seq_rev


col_interval = []
for val in id_nameseq_dict_val:
    col_interval += [val + '_start']
    col_interval += [val + '_end']
interval_df = pd.DataFrame(columns=col_interval)

def variance_calling():
    global pos_set#!
    global seq_seq
    global interval_df
    global title_seq
    bed_columns = ['#contig','start_positioin','end_position','name']
    bed_df = pd.DataFrame(columns=bed_columns)

    for title_seq, seq_seq in single_aln_generator(directory_file_xmfa):
        pos_set = parser_title(title_seq)

        #for interval
        interval_dict = dict()
        for interval,name in zip(pos_set[0],pos_set[1]):
            for val,tag in zip(interval,['_start','_end']):
                interval_dict[name+tag] = int(val)
        interval_df = interval_df.append(interval_dict,ignore_index=True)


        
        if pos_set[1][0]!= REF or int(pos_set[0][0][1]) == 0:#!0 -> 1
            print('WARNING pass aln',pos_set[1][0],pos_set[0][0][1])
            #print('WARNING not contain reference',len(seq_seq[0]),len(seq_seq))
            #print(title_seq)
            continue        

        else:

            #for bed format
            #bed = open(directory_out + '/' +'test'+'.bed','a')
            #bed_df = pd.DataFrame(columns=['#contig','start_positioin','end_position','name'])            
            #bed.write('\t'.join(['#contig','start_positioin','end_position','name'])+'\n')
            bed_dic = {}
            contig,position_real_f = contig_definder(int(pos_set[0][0][0]),find_locus,find_source)
            contig,position_real_s = contig_definder(int(pos_set[0][0][1]),find_locus,find_source)
            for name in pos_set[1]:
                bed_str = [contig] + [str(position_real_f)] + [str(position_real_s)] + [name]

                bed_dic['#contig'] = contig
                bed_dic['start_positioin'] = str(position_real_f)
                bed_dic['end_position'] = str(position_real_s)
                bed_dic['name'] = name
                bed_df = bed_df.append(bed_dic,ignore_index=True)
                #bed.write('\t'.join(bed_str)+'\n')
            #bed.close()           

            if pos_set[2][0] == '-':
                #print('before===',seq_seq)
                seq_seq = seq_reverse(seq_seq)
                #pos_set = list(pos_set)
                #pos_set[0] = [ i[::-1] for i in pos_set[0]]
                #pos_set = tuple(pos_set)
                #print('after===',seq_seq)
                print('Warning reverse strand')                           


            start = int(pos_set[0][0][0])
            end = int(pos_set[0][0][1])
            #with open('log_pos_m.txt','a') as log_pos:
            #    print(file_genome_read[start-1:end] == seq_seq[0].replace('-',''),'\n',
            #          len(file_genome_read[start:end+1]),'\n',
            #          len(seq_seq[0].replace('-','')),'\n',
            #          file=log_pos)
            for dif_set in diffinder(seq_seq,int(pos_set[0][0][0]),seq_seq[0]):
                #print(dif_set)
                #variance
                name_str_dict,variance,position,info = dif_process(dif_set[1],dif_set[0],dif_set[2])
                
                #if pos_set[2][0] == '-':
                    #print(position,int(pos_set[0][0][0])- (int(position)-int(pos_set[0][0][0])))
                    #position = int(pos_set[0][0][0])- (int(position)-int(pos_set[0][0][0]))

                contig,position_real = contig_definder(position,find_locus,find_source)
                #contig = pos_set[1][0]

                with open('alt_variance_dict2.txt','a') as alt_variance_dict_file:
                    print(info,name_str_dict,variance,position,position_real,file=alt_variance_dict_file)
                bin_var = '\t'.join(name_str_dict.values())
                #!position_real
                columns_vcf =[contig,str(position_real),'.',variance[0],','.join(variance[1:]),
                    '40','PASS',info,'GT',bin_var]                  
                
                columns_vcf = '\t'.join(columns_vcf)
                columns_vcf+='\n'
                with open(directory_out +'/'+ name_vcf,'a') as file_vcf:
                        file_vcf.write(columns_vcf)

    interval_df.to_csv(directory_out +'/'+name_vcf[:-4]+'_'+'interval.csv',sep='\t')

    #bed_df.index = bed_df['#contig']
    #del bed_df['#contig']
    bed_df.sort_values(by=['#contig','start_positioin'],ascending=[False,True])
    bed_df.to_csv(directory_out + '/' +name_vcf[:-4]+'.bed',sep='\t',index=False)


variance_calling()
#sorting
if sort:
    with open(directory_out +'/'+name_vcf) as vcf:
        n_comment = 0
        for line in vcf:
            if '##' in line:
               n_comment += 1
    df_vcf = pd.read_csv(directory_out +'/'+name_vcf,sep='\t',header=n_comment)
    df_vcf = df_vcf.sort_values(by=['#CHROM','POS'],ascending=[False,True])

    if delete_ref:
        del df_vcf[REF]
    os.remove(directory_out +'/'+name_vcf)    
    with open(directory_out +'/'+name_vcf,'a') as vcf:
        vcf.write(tag)
    df_vcf.to_csv(directory_out +'/'+name_vcf,sep='\t',index=False,mode='a')

    #df_vcf[df_vcf['INFO']!='WARNING maybe wrong alignment'].to_csv(name_vcf_simple,sep='\t',index=False)

print('end')
#directory_file_xmfa = 
def aln_getter(query_pos,start_inter=100,end_inter=100): 
    contig,position_real = contig_definder(query_pos,find_locus,find_source)
    position_real = query_pos
    for title_seq, seq_seq in single_aln_generator(directory_file_xmfa):
        #print(title_seq)

        pos_set = parser_title(title_seq)

        if pos_set[1][0]!= REF or int(pos_set[0][0][1]) == 0:
            print('WARNING pass aln',pos_set[1][0],pos_set[0][0][1])
            continue
        elif pos_set[1][0]== REF:
            #print(int(pos_set[0][0][0]), query_pos,int(pos_set[0][0][1]))
            if int(pos_set[0][0][0]) <= position_real <= int(pos_set[0][0][1]):
                pos = position_real - int(pos_set[0][0][0])
                pos_gapless = 0
                pos_full = 0
                for i in seq_seq[0]:#range(len(seq_seq[0])):
                    if i != '-':
                        pos_gapless += 1
                    if pos_gapless == pos or pos==0:#0
                        break
                    pos_full += 1 
                start_inter = pos_full - start_inter
                end_inter = pos_full + end_inter
                #print('pos_full=',pos_full,'pos_gapless=',pos_gapless)
                if start_inter<0:
                    start_inter = 0
                if end_inter > len(seq_seq[0])-1:
                    end_inter = len(seq_seq[0])-1
                fast_opened = open(contig + '_' +str(position_real) + '_variance' + '.fna','a')

                for n_seq in range(len(seq_seq)):
                    head = ' '.join(['>',str(position_real),'start_inter=',str(start_inter),'end_inter',str(end_inter),pos_set[1][n_seq],title_seq[n_seq],'\n'])
                    #head = ' '.join(['>',str(position_real),'start_inter=',str(start_inter),'end_inter',str(end_inter),'\n'])

                    fast_opened.write(head)
                    fast_opened.write(seq_seq[n_seq][start_inter:end_inter] + '\n')
                    #print('position_real=',position_real,'query_pos=',query_pos,
                    #    'start_inter=',start_inter,'end_inter=',end_inter,'len(seq_seq[0])-1',len(seq_seq[0])-1,
                    #    pos_set[1][n_seq])
                fast_opened.close()

#aln_getter(164380,start_inter=300,end_inter=300)
