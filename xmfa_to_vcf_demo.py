#By Grigorii Pechkovskii
print('start')

import re
import os

#import linecache
#import tokenize

files = os.listdir()
directory = os.getcwd()
directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/test_mini.xmfa'
#directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/parsnp_edit.xmfa'


file_xmfa = open(directory_file_xmfa)

def get_mauve_index(file_xmfa):
    '''Iter on xmfa header with #
        and return dict with id(key) and name genome(values)
    '''
    id_nameseq_dict = dict()
    for file_line in file_xmfa:
        #FormatVersion Mauve1
        if '#' in file_line:
            find = re.search(r'Sequence(\d*)File\t(.*)\.',file_line)
            if find != None:
                id_nameseq_dict[find.group(1)] = os.path.basename(find.group(2))
        #if '>' in file_line:
        else:
            break
    file_xmfa.close()
    return id_nameseq_dict

id_nameseq_dict = get_mauve_index(file_xmfa)


def head_vcf(name_seq:list):
    '''Function for make vcf header,
       get name_seq: list of name from xmfa header with strict order
    '''
    columns_name =['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    columns_name += name_seq
    head = '\t'.join(columns_name) + '\n'
    with open('test_mini.vcf','w') as file_vcf:
        file_vcf.write(head)
    print(head)

head_vcf(list(id_nameseq_dict.values()))

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
    '''Parsing > line
       get title_seq(> line)
       return positioin start - end sequence
       and name_idseq - name sequece from header name and id in > line
    '''
    position_every = []
    position_every_dict = dict()
    position_every_num_dict = dict()
    name_every = []
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

    name_seq_title = dict()
    for key,val in position_every_num_dict.items() :
        #print(id_nameseq_dict[key])
        name_seq_title[id_nameseq_dict[key]] = val#!id_nameseq_dict frome up namespace
     
    pos = list(name_seq_title.values())
    name_idseq = list(name_seq_title)

    #return name_every,position_every,position_every_dict,position_every_num_dict,pos,name_idseq
    return pos,name_idseq


def diffinder(seq_seq,position,name_seq,head):
    '''Fun fo detect dif betwen aln
    '''
    ref_seq = seq_seq[0]#!reference sequence for compare  
    pos_vcf = int(position[0][0]) #- 1#reference start position
    name_bin_dict = {x: "NA" for x in head}
    

    for sym_num in range(len(ref_seq)):
        
        sym_eq = sym_num
        alt = set()#set alt variance
        alt_dict = dict()
        var_bin = list()
        var_bin_dict = dict()
        alt_num = 0

        if ref_seq[sym_num] != '-':
                #print(ref_seq[sym_num])                
                pos_vcf += 1
        for seq_num,name in zip(range(len(seq_seq)),name_seq):
            #print(name)
            seq = seq_seq[seq_num]  #!Letter, for more clarity code 

            if ref_seq[sym_num] == seq[sym_num]:

                alt_dict[ref_seq[sym_num]] = 0
                var_bin.append(alt_dict[seq[sym_num]])
                #var_bin_dict[name] = alt_dict[seq[sym_num]]
                name_bin_dict[name] = str(alt_dict[seq[sym_num]])

                #var_bin.append(0)                
            else:
                alt.add(seq[sym_num])#set alt variance
                
                if seq[sym_num] not in alt_dict:
                    alt_num += 1
                    alt_dict[seq[sym_num]] = alt_num

                var_bin.append(alt_dict[seq[sym_num]])
                #var_bin_dict[name] = alt_dict[seq[sym_num]]
                name_bin_dict[name] = str(alt_dict[seq[sym_num]])

        if len(alt) == 1:
            sym_eq = sym_num
        else:
            #sym_dif += [sym_num]
            pass
        #under function maybe
        if len(alt) != 0:

            columns_name =['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
            name_bin = '\t'.join(list(name_bin_dict.values()))
            columns_name_dif =['#CHROM',str(pos_vcf),'ID',ref_seq[sym_num],','.join(list(alt_dict)[1:]),'QUAL','FILTER','INFO','FORMAT',name_bin]
            columns_name_dif = '\t'.join(columns_name_dif)
            columns_name_dif += '\n'
            with open('test_mini.vcf','a') as file_vcf:
                file_vcf.write(columns_name_dif)

            yield pos_vcf,ref_seq[sym_num],alt,alt_dict,var_bin,name_bin_dict,columns_name_dif
            #yield 


'''
for i in zip(*seq_seq):
    print(i)
'''
for title_seq, seq_seq in single_aln_generator(directory_file_xmfa):
    #print(title_seq,seq_seq)
    pos_set = parser_title(title_seq)
    #stat = diffinder(seq_seq,pos[1])
    for i in diffinder(seq_seq,pos_set[0],pos_set[1],list(id_nameseq_dict.values())):
        pass

        #print(i)
#req add position parametr if need break loop through equal symbol
sym_seq_start = 'blank'
sym_num_start = 0
sym_seq_end_last = ''
sym_seq_lst = list()
ref_pos = int(pos_set[0][0][0])
for sym_num,sym_seq in enumerate(zip(*seq_seq)):
    sym_num += 1 
    if len(set(sym_seq)) > 1:

        sym_seq_end = sym_seq
        sym_num_end = sym_num

        sym_seq_lst += [list(sym_seq)]
        #print('start=',sym_num_start,sym_seq_start)
        #print('if',sym_num,set(sym_seq),sym_seq_lst)


    elif len(set(sym_seq)) <= 1 and sym_seq_end != sym_seq_end_last:
    #else:#если равны
        
        print('return','start_seq=',sym_seq_start,'\n','start_num',ref_pos + sym_num_start,'\n',
            'end_num=',ref_pos + sym_num_end,'\n','seq_end=',sym_seq_end,'\n','seq_lst=',sym_seq_lst,end='\n')        
        sym_seq_end_last = sym_seq_end

    else:
        sym_num_start = sym_num
        sym_seq_start = sym_seq
        sym_seq_lst = []

def join(massive_str):
    massive_list = []
    for one_str in massive_str:
        massive_list.append(''.join(list(one_str)))
    return massive_list

#seq_seq[3][79:83]
#seq_seq[0][79:83]

print('end')
