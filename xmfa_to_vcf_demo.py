#By Grigorii Pechkovskii
print('start')

import re
import os

#import linecache
#import tokenize

files = os.listdir()
directory = os.getcwd()
directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf_demo/mini_test.xmfa'

#with open(directory_file_xmfa) as file_xmfa:
#    pass



file_xmfa = open(directory_file_xmfa)

def get_mauve_index(file_xmfa):
    number_seq_dict = dict()
    for file_line in file_xmfa:
        #FormatVersion Mauve1
        if '#' in file_line:
            find = re.search(r'Sequence(\d*)File\t(.*)\.',file_line)
            if find != None:
                number_seq_dict[find.group(1)] = os.path.basename(find.group(2))
        #if '>' in file_line:
        else:
            break
    file_xmfa.close()
    return number_seq_dict

number_seq_dict = get_mauve_index(file_xmfa)

seq_dict = dict()



def single_aln_generator(directory_file_xmfa):
    ''' '''
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
                tmp += line.strip() 

        if '=' in line:
            seq_seq.append(tmp)
            tmp = ''
            notfirst = False
            yield title_seq , seq_seq
            title_seq = []
            seq_seq = []

    file_xmfa.close()


def parser_title(title_seq:list):
    ''' '''
    position_every = []
    position_every_dict = dict()
    name_every = []
    for title in title_seq:
        #!!!MAUVE
        name_search = re.search(r'(/.*)\.',title)
        name = name_search.group(1)
        name = os.path.basename(name)
        name_every.append(name)

        position_search = re.search(r'\d*:(\d*)-(\d*)\b',title)
        position = [position_search.group(1),position_search.group(2)]
        position_every += [position]

        position_every_dict[name] = position_every

    return name_every,position_every,position_every_dict


def diffinder(seq_seq,position):
    '''Fun fo detect dif betwen aln
    '''
    ref_seq = seq_seq[0]
    pos_vcf = int(position[0][0]) - 1    

    for sym_num in range(len(ref_seq)):

        alt = set()#множество альтернативных вариантов
        alt_dict = dict()
        var_bin = list()
        alt_num = 0

        if ref_seq[sym_num] != '-':
                #print(ref_seq[sym_num])                
                pos_vcf += 1
        for seq_num in range(len(seq_seq)):

            seq = seq_seq[seq_num]  #! for more clarity code 

            if ref_seq[sym_num] == seq[sym_num]:
                alt_dict[ref_seq[sym_num]] = 0
                var_bin.append(alt_dict[seq[sym_num]])
                #var_bin.append(0)                
            else:
                alt.add(seq[sym_num])#Запишу множество альтернативных вариантов
                
                if seq[sym_num] not in alt_dict:
                    alt_num += 1
                    alt_dict[seq[sym_num]] = alt_num

                var_bin.append(alt_dict[seq[sym_num]])#0 и 1

        if len(alt) != 0:
            yield pos_vcf,ref_seq[sym_num],alt,alt_dict,var_bin   
    #return ref_seq[sym_num],alt,var_bin,alt_dict,pos_vcf



for title_seq, seq_seq in single_aln_generator(directory_file_xmfa):
    #print(title_seq,seq_seq)
    pos = parser_title(title_seq)
    #stat = diffinder(seq_seq,pos[1])
    for i in diffinder(seq_seq,pos[1]):
        print(i)

print('end')
