print('START PIPLINE')
import os
import subprocess
import re

import numpy as np
import pandas as pd

import pipeline_base

#config
mauve = '/home/strain4/Desktop/content/bioinf_prog/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve'
xmfa_to_vcf = './xmfa_to_vcf_demo.py'
bcftools = '/home/strain4/Desktop/content/bioinf_prog/bcftools/bcftools'
bgzip = 'bgzip'
vcf_merger = './vcf_merger.py'


work_dir = '/home/strain4/Desktop/piplines/genomics_pipline_supply/genome/'
REF = '/home/strain4/Desktop/piplines/genomics_pipline_supply/genome/GCF_000008445.1_ASM844v1_genomic.fna'
name_exp = 'exp_test_B1'
out_dir = '/home/strain4/Desktop/piplines/genomics_pipline_supply/' + name_exp + '/' # / impotant
file_gbk = '/home/strain4/Desktop/piplines/genomics_pipline_supply/' + 'AmesAncestor_GCF_000008445.1.gbk'

#header = 15
BED = True
pipeline_version = 'pipeline version 0.05\n'


os.mkdir(out_dir)

logfile_mauve_path = out_dir+'log_mauve_file.txt'
logfile_mauve = open(logfile_mauve_path,'a')

logfile_path = out_dir+'log_file.txt'
logfile = open(logfile_path,'a')
logfile.write(pipeline_version)

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



find_locus, find_source ,find_source_real = pipeline_base.contig_finder_gbk(file_gbk)
for vcf_path in vcf_path_lst:
    logfile = open(logfile_path,'a')
    logfile.write('\n'+vcf_path+'\n')
    logfile.close()
    logfile = open(logfile_path,'a')
    print('\n'+vcf_path+'\n')
    
    header,vcf_head,type_head = pipeline_base.vcf_head_process(vcf_path)


    vcf_editer_path = vcf_path[0:-4]+'_editer.vcf'
    vcf_editer_opened = open(vcf_editer_path,'a')
    
    vcf = pd.read_csv(vcf_path,sep='\t',header = header)
    vcf.index = vcf['POS']
    #find_locus, find_source ,find_source_real = pipeline_base.contig_finder_gbk(file_gbk)#!
    vcf,vcf_pos_old_new = pipeline_base.position_editer(vcf.copy(),find_locus,find_source,old_new=True)
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

    
    vcf = pd.read_csv(vcf_path,sep='\t',header = header)
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

    #vcf[vcf['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X')]

    pos_will_del_old += dup_list



    logfile = open(logfile_path,'a')
    logfile.write('pos_will_del_old ' + str(pos_will_del_old))#!open
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
#!!!!for 1 do not work
merged_vcf_path = vcf_out + 'merged_' + name_exp + '.vcf'
with open(merged_vcf_path,'a') as file_merged:
    #bcftools_run = subprocess.Popen([bcftools,'merge','--merge','all','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged,stderr=logfile)#,stderr=logfile)
    bcftools_run = subprocess.Popen([bcftools,'merge','--merge','none','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged,stderr=logfile)#,stderr=logfile)

    bcftools_run.wait()


header,vcf_head,type_head = pipeline_base.vcf_head_process(merged_vcf_path)


print('START vcf_merger')
merged_final = vcf_out +  'merged_final_'  + name_exp + '.vcf'
vcf_merger_run = subprocess.Popen([vcf_merger,
        "-v",  merged_vcf_path,
        "-r", REF,        
        "-g", file_gbk,
        "-d", vcf_out,
        "-l", logfile_path,
        "-t", str(header),
        "-i",','.join(intervals_path),
        "-o", merged_final],universal_newlines=True)
vcf_merger_run.wait()


logfile.close()
print('END PIPLINE')