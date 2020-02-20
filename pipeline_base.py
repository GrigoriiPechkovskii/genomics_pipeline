print('start_base')

import os
import subprocess
import re

import numpy as np
import pandas as pd


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
                    #print(line)
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
                            fasta_out_str +=  ' '.join([fasta_name,str(start),str(end)]) + '\n' +  fasta_seq[start:end] + '\n'
                            break
                        else:
                            #print('WARNING sequence from several contig')
                            fasta_out_str += ' '.join([fasta_name,str(start),str(len(fasta_seq))]) + '\n' +  fasta_seq[start:end] + '\n'
                            end-=len(fasta_seq)
                            start= 0
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

            elif sym=='W':
                seq_rev += 'W'
            elif sym=='S':
                seq_rev += 'S'
            elif sym=='M':
                seq_rev += 'K'
            elif sym=='K':
                seq_rev += 'M'
            elif sym=='R':
                seq_rev += 'Y'
            elif sym=='Y':
                seq_rev += 'R'
            elif sym=='B':
                seq_rev += 'V'
            elif sym=='D':
                seq_rev += 'H'
            elif sym=='H':
                seq_rev += 'D'
            elif sym=='V':
                seq_rev += 'B'
            elif sym=='N':
                seq_rev += 'N'
            elif sym=='Z':
                seq_rev += 'Z'
            else:
                print('Error with nuc sym')


        seq_rev = seq_rev[::-1]
        seq_seq_rev += [seq_rev]
    return seq_seq_rev


def contig_finder_gbk(file_gbk_dir):
    ''' '''
    with open(file_gbk_dir) as file_gbk_opened:
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



def position_editer(vcf,find_locus,find_source,old_new=False):
    vcf_pos_old_new = {}
    pos_lst = []
    for position in vcf['POS']:
        contig,position_real = contig_definder(position,find_locus,find_source)
        pos_lst.append(position_real)
        #vcf.loc[position,'POS'] = position_real
        if old_new:
            vcf_pos_old_new[contig+'_'+str(position_real)] = position

    vcf['POS'] = np.array(pos_lst)
    vcf.index = vcf['POS']

    if old_new:
        return vcf,vcf_pos_old_new
    else:
        return vcf


def vcf_head_process(vcf_dir):
    header = 0
    vcf_head = ''
    vcf_opened = open(vcf_dir,'r')
    for line in vcf_opened:
        if line.startswith('##'):
            header += 1
            vcf_head += line
            continue
        elif line.startswith('#'):
            head_vcf_line = line.strip().split('\t')
        else:
            break
    type_head = {i:str for i in head_vcf_line[9:]}

    return header,vcf_head,type_head



def position_editer111(vcf):
    pos_lst = []
    for position in vcf['POS']:
        contig,position_real = contig_definder(position,find_locus,find_source)
        pos_lst.append(position_real)
        #vcf.loc[position,'POS'] = position_real
    vcf['POS'] = np.array(pos_lst)
    vcf.index = vcf['POS']
    return vcf


print('end_base')