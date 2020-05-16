print('START PIPLINE')
import os
import subprocess
import re

import numpy as np
import pandas as pd

import time

import pipeline_base

import config



def mauve_process(genome_grouped,start_num=0,tag=''):   

    #run mauve
    print('mauve_process')
    log_time = []
    xmfa_dir_lst = []
    #vcf_path_lst = []
    for num,genome_quary_list in enumerate(genome_grouped,start_num):

        genome_quary_list_name = [os.path.basename(k) for k in [os.path.splitext(i)[0] for i in genome_quary_list]]

        genome_quary = '|'.join(genome_quary_list_name)

        #name_out = name_exp + '_' + 'group_' + str(num)
        name_out = name_exp + '_' + genome_quary + '_group_' + str(num)

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


        xmfa_dir_lst += [group_file]

        current_time = pipeline_base.timecheck(name_out)
        log_time.append(current_time)
        print(current_time)

    logfile.write(''.join(log_time))

    return xmfa_dir_lst


'''        xmfa_dir = group_file

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

    return vcf_path_lst'''

def xmfa_to_vcf_process(xmfa_dir_lst):

    vcf_path_lst = []
    log_time = []

    for xmfa_dir in xmfa_dir_lst:

        name_out = os.path.basename(xmfa_dir)


        vcf_path = os.path.join(vcf_out + name_out + '.vcf')

        if os.access(vcf_path,os.F_OK):
            print(name_out + '.vcf','exist, pass vcf')
            continue

        else:
            vcf_path_lst += [vcf_path]

            xmfa_to_vcf_run = subprocess.Popen([xmfa_to_vcf,
                "-x",  xmfa_dir,
                "-r", os.path.basename(REF)[:-4],
                "-o", vcf_out,
                "-g", file_gbk,
                "-n", name_out + '.vcf',
                ],universal_newlines=True,stdout=logfile)

            #os.path.dirname(xmfa_dir)
            #vcf_path_lst += [vcf_out + name_out + '.vcf']

            xmfa_to_vcf_run.wait()


            log_time.append(pipeline_base.timecheck(name_out))


    logfile.write(''.join(log_time))


    return vcf_path_lst


def vcf_editer(vcf_path_lst):

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



#def bcftools_merge(vcf_path_lst):
def bgzip_bcftools_indexing(vcf_path_lst):

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
            #bgzip_run = subprocess.Popen([bgzip,vcf_path],universal_newlines=True)

            with open(vcf_path + '.gz','w') as vcf_path_gz_opened:
                bgzip_run = subprocess.Popen([bgzip,'-c',vcf_path],universal_newlines=True,stdout=vcf_path_gz_opened)
                bgzip_run.wait()

            bcftools_index_run = subprocess.Popen([bcftools,'index',vcf_path + '.gz'],universal_newlines=True)
            bcftools_index_run.wait()

    return vcf_path_lst_gz


    #del vcf_path_lst_gz[-1]# without add all
    #err = open("/home/strain4/Desktop/piplene_mauve/work_dir/out/err.txt",'a' )

    #del whole path interval
def intervals_path_finder(vcf_out):
    if BED:
        intervals_path = [vcf_out+i for i in os.listdir(vcf_out) if '.bed' in i]#!! move
    else:
        intervals_path = [vcf_out+i for i in os.listdir(vcf_out) if 'interval' in i]
    return intervals_path

    '''
    for n_intervals in  range(len(intervals_path)):
        if name_exp + '_' + 'group_' + str(num) in intervals_path[n_intervals]:
            del intervals_path[n_intervals]
            break'''

def bcftools_merge(vcf_path_lst_gz):
    #del intervals_path[]
    #!!!!for 1 do not work
    merged_vcf_path = vcf_out + 'merged_' + name_exp + '.vcf'
    with open(merged_vcf_path,'a') as file_merged:
        #bcftools_run = subprocess.Popen([bcftools,'merge','--merge','all','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged,stderr=logfile)#,stderr=logfile)
        bcftools_run = subprocess.Popen([bcftools,'merge','--merge','none','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged,stderr=logfile)#,stderr=logfile)

        bcftools_run.wait()

    return merged_vcf_path

#intervals_path, merged_vcf_path = bcftools_merge() #intervals_path not here!!!



def vcf_merger_run(intervals_path, merged_vcf_path):
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





#config

mauve = config.mauve
xmfa_to_vcf = config.xmfa_to_vcf
bcftools = config.bcftools
bgzip = config.bgzip
vcf_merger = config.vcf_merger

work_dir = config.work_dir  #genome dir
REF = config.REF
name_exp = config.name_exp
out_dir = config.out_dir    #aln dir
file_gbk = config.file_gbk

files_add = config.files_add



adding_mode = False
mauve_aln_exist = False
xmfa_to_vcf_process_exist = False


if adding_mode:
    work_dir = config.files_add
    REF = config.REF_add


#header = 15
BED = True
pipeline_version = 'pipeline version 0.05\n'


if not os.access(out_dir,os.F_OK):
    os.mkdir(out_dir)


vcf_out = out_dir + 'vcf_out/'#os.path.join(out_dir,'vcf_out')!!!!!
if not os.access(vcf_out,os.F_OK):
    os.mkdir(vcf_out)

#logfile_mauve_path = out_dir+'log_mauve_file.txt'
logfile_mauve_path = os.path.join(out_dir,'log_mauve_file.txt')

logfile_mauve = open(logfile_mauve_path,'a')

#logfile_path = out_dir+'log_file.txt'
logfile_path = os.path.join(out_dir,'log_file.txt')

logfile = open(logfile_path,'a')
logfile.write(pipeline_version)


def grouping_genome(work_dir,REF):

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

    return genome_grouped

#genome_grouped += [[REF] + files_genome]

#del genome_grouped[-1]#!
#del genome_grouped[-1]#!


'''
if mauve_aln_exist:
    xmfa_dir_lst = [os.path.join(out_dir,dir_group,dir_group) for dir_group in os.listdir(out_dir) if 'group' in dir_group]


vcf_path_lst_exist = []
if xmfa_to_vcf_process_exist:
    #vcf_path_lst = [xmfa_dir + '.vcf' for xmfa_dir in xmfa_dir_lst]
    vcf_path_lst_exist = [os.path.join(vcf_out, dir_group + '.vcf') for dir_group in os.listdir(out_dir) if 'group' in dir_group]



xmfa_dir_lst = mauve_process()

vcf_path_lst = xmfa_to_vcf_process(xmfa_dir_lst)
vcf_editer(vcf_path_lst)

intervals_path, merged_vcf_path = bcftools_merge(vcf_path_lst_exist + vcf_path_lst)
#intervals_path, merged_vcf_path = bcftools_merge(vcf_path_lst)

vcf_merger_run()
'''


def pipeline_full():
    adding_mode = False
    mauve_aln_exist = False
    xmfa_to_vcf_process_exist = False
    vcf_path_lst_exist = []

    genome_grouped = grouping_genome(work_dir,REF)

    xmfa_dir_lst = mauve_process(genome_grouped)

    vcf_path_lst = xmfa_to_vcf_process(xmfa_dir_lst)
    vcf_editer(vcf_path_lst)

    vcf_path_lst_gz = bgzip_bcftools_indexing(vcf_path_lst)
    intervals_path = intervals_path_finder(vcf_out)

    merged_vcf_path = bcftools_merge(vcf_path_lst_exist + vcf_path_lst_gz)

    vcf_merger_run(intervals_path, merged_vcf_path)
 
def pipeline_part_aln():
    adding_mode = False
    mauve_aln_exist = False
    xmfa_to_vcf_process_exist = False

    genome_grouped = grouping_genome(work_dir,REF)

    xmfa_dir_lst = mauve_process(genome_grouped)
   

def pipeline_adding_aln():

    adding_mode = True

    if adding_mode:
        work_dir = config.files_add
        REF = config.REF_add

    genome_grouped = grouping_genome(work_dir,REF)

    xmfa_dir_lst = mauve_process(genome_grouped)
    

def pipelene_adding_vcf():
    mauve_aln_exist = True

    if mauve_aln_exist:
        xmfa_dir_lst = [os.path.join(out_dir,dir_group,dir_group) for dir_group in os.listdir(out_dir) if 'group' in dir_group]
    
    vcf_path_lst = xmfa_to_vcf_process(xmfa_dir_lst)
    vcf_editer(vcf_path_lst)

    vcf_path_lst_gz = bgzip_bcftools_indexing(vcf_path_lst)


def pipelene_vcf_merger():

    intervals_path = intervals_path_finder(vcf_out)
    vcf_path_lst_gz = [os.path.join(vcf_out, dir_group + '.vcf.gz') for dir_group in os.listdir(out_dir) if 'group' in dir_group]

    merged_vcf_path = bcftools_merge(vcf_path_lst_gz)

    vcf_merger_run(intervals_path, merged_vcf_path)

def pipeline_adding_aln_vcf():
    pipeline_adding_aln()
    pipelene_adding_vcf()
    pass

time_block_start = time.time()
time_initial = time.time()

#pipeline_full()
#pipeline_adding_aln()
#pipelene_adding_vcf()
#pipelene_vcf_merger()

pipeline_part_aln()

logfile.close()
print('END PIPLINE')


#time!!