#By Grigorii Pechkovskii
'''
Check ext with blank alt
Check positioins shift 1
'''
print('start')

import re
import os

#import linecache
#import tokenize

files = os.listdir()
directory = os.getcwd()
directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/test_mini.xmfa'
directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/parsnp_edit.xmfa'
directory_file_xmfa = '/home/strain4/Desktop/xmfa_to_vcf/mauve_out1.xmfa'
name_vcf = 'test_mini2.vcf'

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
id_nameseq_dict_val = list(id_nameseq_dict.values())

def head_vcf(name_seq:list):
    '''Function for make vcf header,
       get name_seq: list of name from xmfa header with strict order
    '''
    columns_name =['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    columns_name += name_seq
    head = '\t'.join(columns_name) + '\n'
    with open(name_vcf,'w') as file_vcf:
        file_vcf.write(head)
    print(head)

head_vcf(id_nameseq_dict_val)

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


def diffinder(seq_seq):
    #req add position parametr if need break loop through equal symbol
    sym_seq_start = 'blank'
    sym_seq_end = ''
    sym_num_start = 0
    sym_seq_lst = list()
    ref_pos = int(pos_set[0][0][0])
    ref_seq = seq_seq[0]#!reference sequence for compare  
    pos_vcf = -1
    start_pos = ref_pos
    first_flag = True
    for sym_num,sym_seq in enumerate(zip(*seq_seq)):

        if ref_seq[sym_num] != '-':
                pos_vcf += 1

        if len(set(sym_seq)) > 1:#input in loop if alternative exist
            sym_seq_end = sym_seq
            sym_num_end = sym_num
            sym_seq_lst += [list(sym_seq)]

            if first_flag:
                start_pos = ref_pos + pos_vcf
                first_flag = False

        elif len(set(sym_seq)) == 1 and len(sym_seq_lst) != 0:#here break
            sym_seq_lst = [list(sym_seq_start)] + sym_seq_lst

            yield start_pos,sym_seq_lst

            first_flag = True
            sym_num_start = sym_num
            sym_seq_start = sym_seq
            sym_seq_lst = []

        else:
            sym_seq_start = sym_seq
        
def join(massive_str):
    massive_list = []
    for one_str in massive_str:
        massive_list.append(''.join(list(one_str)))
    return massive_list

def dif_process(seq_lst,position):
    variance_bin_dict = dict()

    name_str_dict = {x: "NA" for x in id_nameseq_dict_val}
    variance_lst = join(list((zip(*seq_lst))))
    variance_set = set((zip(*seq_lst)))

    alt_variance_set  = set()
    alt_variance_dict = dict()
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
        if name_str_dict[name] != 'NA':
            name_str_dict[name] = alt_variance_dict[name_str_dict[name]]
    
    #for SNP
    alt_variance_dict_n = {}

    if all(len(s) == 2 for s in list(alt_variance_dict.keys())) and all('-' not in s for s in list(alt_variance_dict.keys())):        
        for alt_key,alt_val in alt_variance_dict.items():
            alt_variance_dict_n[alt_key[1:]] = alt_val
        position += 1
        return name_str_dict,alt_variance_dict_n,position
    
    dash_logic = True
    if dash_logic:
        for alt_key,alt_val in alt_variance_dict.items():
            alt_variance_dict_n[alt_key.replace('-','')] = alt_val
        return name_str_dict,alt_variance_dict_n,position

    return name_str_dict,alt_variance_dict,position



file_gbk = directory + '/AmesAncestor_GCF_000008445.1.gbk'

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
            position_real = position - source[0]+ 1
            return locus,position_real

for title_seq, seq_seq in single_aln_generator(directory_file_xmfa):
    pos_set = parser_title(title_seq)
    if pos_set[1][0]!= 'AmesAncestor_GCF_000008445.1':
        print('if')
        continue
        
    else:
        #print('else')
        for dif_set in diffinder(seq_seq):
            name_str_dict,alt_variance_dict,position = dif_process(dif_set[1],dif_set[0])

            bin_var = '\t'.join(name_str_dict.values())

            contig,position_real = contig_definder(position,find_locus,find_source)

            columns_vcf =[contig,str(position),'.',list(alt_variance_dict)[0],','.join(list(alt_variance_dict)[1:]),
                '40','PASS','NA','GT',bin_var]
            columns_vcf = '\t'.join(columns_vcf)
            columns_vcf+='\n'
            with open(name_vcf,'a') as file_vcf:
                    file_vcf.write(columns_vcf)




#position = 1000
#q = contig_definder(position,find_locus,find_source)

print('end')

