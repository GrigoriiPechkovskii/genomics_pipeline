#!/usr/bin/python3

#By Grigorii Pechkovskii
'''

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
parser.add_argument('-r', '--ref',action='store', help='Reference fasta')
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

    directory_file_xmfa = directory + '/' + 'test_mini.xmfa'
    directory_file_xmfa = '/home/strain4/Desktop/fin_script/xmfa_to_vcf/exp_A2_group_0'

    REF = 'AmesAncestor_GCF_000008445.1'#test_mini
    #REF ='AmesAncestor_GCF_0000084451'
    REF = 'GCF_000008445.1_ASM844v1_genomic'#GI_AAJ_out1
    #REF = 'Ames_Ancestor_ref_GCF_000008445.1_ASM844v1_genomic.fna'#parsnp.xmfa
    directory_out = directory
    index_type = 'mauve'
    #name_vcf_simple = 'test_sim.vcf'
    file_gbk = directory + '/AmesAncestor_GCF_000008445.1.gbk'
    name_vcf = 'test_exp_1_group_1.vcf'

#some important options
sort = True
delete_ref = True
NORM = True

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
                    head = ' '.join(['>' + pos_set[1][n_seq],str(position_real),'start_inter=' + str(start_inter),'end_inter=' + str(end_inter),title_seq[n_seq].replace('>','').replace(' ','_'),'\n'])
                    #head = ' '.join(['>',str(position_real),'start_inter=',str(start_inter),'end_inter',str(end_inter),'\n'])

                    fast_opened.write(head)
                    fast_opened.write(seq_seq[n_seq][start_inter:end_inter] + '\n')
                    #print('position_real=',position_real,'query_pos=',query_pos,
                    #    'start_inter=',start_inter,'end_inter=',end_inter,'len(seq_seq[0])-1',len(seq_seq[0])-1,
                    #    pos_set[1][n_seq])
                fast_opened.close()

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

    def fasta_getter(self,start,end,contig='whole'):
        if contig == 'whole':
            length_seq = 0
            fasta_out_str = ''
            start_out = start
            end_out = end
            start-=1                
            if contig == 'whole':
                for fasta_seq,fasta_name in zip(self.seq_lst,self.name_lst):
                    0<0
                    if start<=len(fasta_seq):# and end<=len(fasta_seq):
                        if end<=len(fasta_seq):
                            #print(fasta_name,start,end)
                            #print(fasta_seq[start:end])
                            fasta_out_str +=  ' '.join([fasta_name,str(start),str(end)]) + '\n' +  fasta_seq[start:end] + '\n'
                            break
                        else:
                            #print('WARNING sequence from several contig')
                            #print(fasta_name,start,len(fasta_seq))
                            #print(fasta_seq[start:end])
                            fasta_out_str += ' '.join([fasta_name,str(start),str(len(fasta_seq))]) + '\n' +  fasta_seq[start:end] + '\n'
                            end-=len(fasta_seq)
                            start= 0
                            #print('end=', end)
                    else:
                        start-= len(fasta_seq)
                        end-=len(fasta_seq)
                        pass
            else:
                for fasta_seq,fasta_name in zip(self.seq_lst,self.name_lst):
                    if contig in fasta_name:
                        fasta_seq[start,end]


        fasta_out = open(contig +' '+ str(start_out) +':'+ str(end_out) + '.fna','w')
        fasta_out.write(fasta_out_str)
        fasta_out.close()

        return fasta_out_str

file_fasta = 'test_mini_vcf.fna'

fasta = SequenceFasta(file_fasta)
fasta.seq_process(strip=True)
#sequence = ''.join(fasta.seq_lst)


def fasta_getter(fasta_path,start=0,end=1):
    fasta_opened = open(fasta_path)
    for line in fasta_opened:
        if line.startswith('>'):
            print(line)
        pass

#aln_getter(164380,start_inter=300,end_inter=300)

print('end')
