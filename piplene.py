print('START PIPLINE')
import os
import subprocess
import re

import numpy as np
import pandas as pd

#config
mauve = '/home/strain4/Desktop/content/bioinf_prog/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve'
xmfa_to_vcf = '/home/strain4/Desktop/fin_script/xmfa_to_vcf/xmfa_to_vcf_demo.py'
bcftools = '/home/strain4/Desktop/content/bioinf_prog/bcftools/bcftools'
bgzip = 'bgzip'
vcf_merger = '/home/strain4/Desktop/fin_script/genomics_pipline/vcf_merger.py'

'''
work_dir = '/home/strain4/Desktop/piplene_mauve/work_dir/'
REF = '/home/strain4/Desktop/piplene_mauve/work_dir/x.fna'
name_exp = 'test'
out_dir = "/home/strain4/Desktop/piplene_mauve/out_test4/"
'''

work_dir = '/home/strain4/Desktop/fin_script/test_genomics_pipline/genome/'
REF = '/home/strain4/Desktop/fin_script/test_genomics_pipline/genome/GCF_000008445.1_ASM844v1_genomic.fna'
name_exp = 'exp_A3'
out_dir = '/home/strain4/Desktop/fin_script/test_genomics_pipline/expA3_test/' # / impotant
file_gbk = '/home/strain4/Desktop/fin_script/test_genomics_pipline/AmesAncestor_GCF_000008445.1.gbk'
header = 15
BED = True


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
    vcf_pos_old_new = {}
    pos_lst = []
    for position in vcf['POS']:
        contig,position_real = contig_definder(position,find_locus,find_source)
        pos_lst.append(position_real)
        #vcf.loc[position,'POS'] = position_real
        vcf_pos_old_new[contig+'_'+str(position_real)] = position

    vcf['POS'] = np.array(pos_lst)
    vcf.index = vcf['POS']
    return vcf,vcf_pos_old_new


os.mkdir(out_dir)

logfile_mauve_path = out_dir+'log_mauve_file.txt'
logfile_mauve = open(logfile_mauve_path,'a')

logfile_path = out_dir+'log_file.txt'
logfile = open(logfile_path,'a')

files = os.listdir(work_dir)
files_genome = [work_dir + file for file in files if file.endswith(('.fa','.fasta','.fna'))]

files_genome.remove(REF)
#for grouping
amount = 1
count = 0
genome_run = []
genome_grouped = []
for genome in files_genome:
    genome_run += [genome]
    count +=1        
    if amount == count:
        genome_grouped += [[REF] + genome_run]
        genome_run = []
        count = 0

#genome_grouped += [[REF] + files_genome]

#del genome_grouped[-1]#!
#del genome_grouped[-1]#!

vcf_out = out_dir + 'vcf_out/'
os.mkdir(vcf_out)

#run mauve
print('run mauve')
vcf_path_lst = []
for num,genome_quary_list in enumerate(genome_grouped):
    genome_quary = ' '.join(genome_quary_list)

    name_out = name_exp + '_' + 'group_' + str(num)

    logfile = open(logfile_path,'a')
    logfile.write('\n'+name_out+'\n')#!
    logfile.close()
    logfile = open(logfile_path,'a')
    print('\n'+name_out+'\n')

    group_dir =  out_dir + name_out
    group_file = group_dir + '/' + name_out
    os.mkdir(group_dir)

    mauve_run = subprocess.Popen([mauve,"--output=" + group_file, *genome_quary_list],universal_newlines=True,stdout=logfile_mauve)
    mauve_run.wait()

    xmfa_dir = group_file

    xmfa_to_vcf_run = subprocess.Popen([xmfa_to_vcf,
        "-x",  xmfa_dir,
        "-r", os.path.basename(REF)[:-4],
        "-o", vcf_out,
        "-g", file_gbk,
        "-n", name_out + '.vcf',
        ],universal_newlines=True,stdout=logfile)

    os.path.dirname(group_file)
    vcf_path_lst += [vcf_out + name_out + '.vcf']

    xmfa_to_vcf_run.wait()



find_locus, find_source ,find_source_real = contig_finder_gbk(file_gbk)
for vcf_path in vcf_path_lst:
    logfile = open(logfile_path,'a')
    logfile.write('\n'+vcf_path+'\n')
    logfile.close()
    logfile = open(logfile_path,'a')
    print('\n'+vcf_path+'\n')
    
    header2 = 0
    vcf_head = ''
    vcf_opened = open(vcf_path)
    for line in vcf_opened:
        if '##' in line:
            vcf_head += line
            header2 += 1
    vcf_opened.close()

    vcf_editer_path = vcf_path[0:-4]+'_editer.vcf'
    vcf_editer_opened = open(vcf_editer_path,'a')
    
    vcf = pd.read_csv(vcf_path,sep='\t',header = header2)
    vcf.index = vcf['POS']
    #find_locus, find_source ,find_source_real = contig_finder_gbk(file_gbk)#!
    vcf,vcf_pos_old_new = position_editer(vcf.copy())
    vcf_editer_opened.write(vcf_head)    
    vcf_editer_opened.close()
    vcf.to_csv(vcf_editer_path,sep='\t',index=False,mode='a')

    vcf_err_path = vcf_path[:-4]+'_err.vcf'
    vcf_out_path = vcf_path[:-4]+'_out.vcf'
    vcf_err_opened  = open(vcf_err_path,'a')
    vcf_out_opened = open(vcf_out_path,'a')
    bcftools_norm_run = subprocess.Popen([bcftools,'norm','--check-ref', 'w', '-f',REF,vcf_editer_path],stdout=vcf_out_opened,stderr=vcf_err_opened)
    bcftools_norm_run.wait()
    vcf_err_opened.close()
    vcf_out_opened.close()
    
    vcf_err_df = pd.read_csv(vcf_err_path,sep='\t',names=['type','contig','POS','ref','alt'])

    vcf_err_df.drop(vcf_err_df[vcf_err_df['contig'].isna()].index,inplace=True)
    pos_will_del =(vcf_err_df[0:-1]['contig']  +'_'+ vcf_err_df[0:-1]['POS'].astype(int).astype(str)).values
    
    pos_will_del_old = []
    for pos_del in pos_will_del:
        pos_will_del_old.append(vcf_pos_old_new[pos_del])

    
    vcf = pd.read_csv(vcf_path,sep='\t',header = header2)
    vcf.index = vcf['POS']    

    logfile = open(logfile_path,'a')
    logfile.write('\n'+str(vcf[vcf.duplicated(['POS'],keep=False)])+'\n')
    logfile.close()
    logfile = open(logfile_path,'a')

    print('\n'+str(vcf[vcf.duplicated(['POS'],keep=False)])+'\n',file=logfile)
    dup_list = list(vcf[vcf.duplicated(['POS'],keep=False)]['POS'].values)
    #logfile.write(str(pos_will_del_old))

    vcf.drop(pos_will_del_old,inplace=True)    
    vcf.drop_duplicates(['POS'],keep=False, inplace=True)

    pos_will_del_old += dup_list

    logfile = open(logfile_path,'a')
    logfile.write(str(pos_will_del_old))#!open
    logfile.close()
    logfile = open(logfile_path,'a')
    

    vcf_opened = open(vcf_path,'w')
    vcf_opened.write(vcf_head)    
    vcf_opened.close()
    vcf.to_csv(vcf_path,sep='\t',index=False,mode='a')






vcf_path_lst_gz = [vcf_path + '.gz' for vcf_path in vcf_path_lst]

with_norm = False
if with_norm:
    for vcf_gz,vcf_path in zip(vcf_path_lst_gz,vcf_path_lst):

        #bgzip_run = subprocess.Popen([bgzip,vcf_path],universal_newlines=True)
        #bgzip_run.wait()

        #bcftools_index_run = subprocess.Popen([bcftools,'index',vcf_gz],universal_newlines=True)
        #bcftools_index_run.wait()

        #work without gzip

        bcftools_norm_run = subprocess.Popen([bcftools,'norm','-Oz','--check-ref', 'xw', '-f',REF,vcf_path,'-o',vcf_gz],universal_newlines=True)#,stderr=logfile)
        #bcftools_norm_run = subprocess.Popen([bcftools,'norm','-Oz','-f',REF,vcf_gz,'-o',vcf_gz],universal_newlines=True)#,stderr=logfile)
        bcftools_norm_run.wait()

        bcftools_index_run = subprocess.Popen([bcftools,'index','-f',vcf_gz],universal_newlines=True)
        bcftools_index_run.wait()
elif not with_norm:
    for vcf_path in vcf_path_lst:
        bgzip_run = subprocess.Popen([bgzip,vcf_path],universal_newlines=True)
        bgzip_run.wait()

        bcftools_index_run = subprocess.Popen([bcftools,'index',vcf_path + '.gz'],universal_newlines=True)
        bcftools_index_run.wait()

#del vcf_path_lst_gz[-1]# without add all
#err = open("/home/strain4/Desktop/piplene_mauve/work_dir/out/err.txt",'a' )

#del whole path interval

if BED:
    intervals_path = [vcf_out+i for i in os.listdir(vcf_out) if '.bed' in i]
else:
    intervals_path = [vcf_out+i for i in os.listdir(vcf_out) if 'interval' in i]

'''
for n_intervals in  range(len(intervals_path)):
    if name_exp + '_' + 'group_' + str(num) in intervals_path[n_intervals]:
        del intervals_path[n_intervals]
        break'''

#del intervals_path[]

merged_vcf_path = vcf_out + 'merged.vcf'
with open(merged_vcf_path,'a') as file_merged:
    bcftools_run = subprocess.Popen([bcftools,'merge','--merge','all','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged)#,stderr=logfile)
    bcftools_run.wait()



print('START vcf_merger')
merged_final = vcf_out +  'merged_final.vcf'
vcf_merger_run = subprocess.Popen([vcf_merger,
        "-v",  merged_vcf_path,
        "-r", REF,        
        "-g", file_gbk,
        "-d", vcf_out,
        "-l", logfile_path,
        "-t", '15',
        "-i",','.join(intervals_path),
        "-o", merged_final        
        ],universal_newlines=True)
vcf_merger_run.wait()

'''parser.add_argument('-v', '--vcf',action='store', help='File vcf')
parser.add_argument('-r', '--ref',action='store', help='Reference fasta')
parser.add_argument('-g', '--gbk-file',action='store', help='File gbk')
parser.add_argument('-d', '--dir',action='store', help='Work directory')
parser.add_argument('-l', '--log',action='store', help='Work directory')
parser.add_argument('-h', '--header-vcf',action='store',default=15, help='Header vcf')'''

logfile.close()
print('END PIPLINE')