print('START PIPLINE')
import os
import subprocess

#config
mauve = '/home/strain4/Desktop/content/bioinf_prog/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve'
xmfa_to_vcf = '/home/strain4/Desktop/xmfa_to_vcf/xmfa_to_vcf_demo.py'
bcftools = "/home/strain4/Desktop/content/bioinf_prog/bcftools/bcftools"
bgzip = 'bgzip'

'''
work_dir = '/home/strain4/Desktop/piplene_mauve/work_dir/'
REF = '/home/strain4/Desktop/piplene_mauve/work_dir/x.fna'
name_exp = 'test'
out_dir = "/home/strain4/Desktop/piplene_mauve/out_test4/"
'''

work_dir = '/home/strain4/Desktop/piplene_mauve/exp1/'
REF = '/home/strain4/Desktop/piplene_mauve/exp1/GCF_000008445.1_ASM844v1_genomic.fna'
name_exp = 'exp_1'
out_dir = "/home/strain4/Desktop/piplene_mauve/exp_test_for4_with_norm2/" # / impotant


os.mkdir(out_dir)

logfile_name = work_dir+'logfile.txt'
logfile = open(logfile_name,'a')

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

genome_grouped += [[REF] + files_genome]


vcf_out = out_dir + 'vcf_out/'
os.mkdir(vcf_out)

#run mauve
vcf_path_lst = []
for num,genome_quary_list in enumerate(genome_grouped):
    genome_quary = ' '.join(genome_quary_list)

    name_out = name_exp + '_' + 'group_' + str(num)

    group_dir =  out_dir + name_out
    group_file = group_dir + '/' + name_out
    os.mkdir(group_dir)

    mauve_run = subprocess.Popen([mauve,"--output=" + group_file, *genome_quary_list],universal_newlines=True,stdout=logfile)
    mauve_run.wait()

    xmfa_dir = group_file

    xmfa_to_vcf_run = subprocess.Popen([xmfa_to_vcf,
        "-x",  xmfa_dir,
        "-r", os.path.basename(REF)[:-4],
        "-o", vcf_out,
        "-g", "/home/strain4/Desktop/xmfa_to_vcf/AmesAncestor_GCF_000008445.1.gbk",
        "-n", name_out + '.vcf',
        ],universal_newlines=True,stdout=logfile)

    os.path.dirname(group_file)
    vcf_path_lst += [vcf_out + name_out + '.vcf']

    xmfa_to_vcf_run.wait()

#del vcf_path_lst[-1]# without add all

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

del vcf_path_lst_gz[-1]# without add all
#err = open("/home/strain4/Desktop/piplene_mauve/work_dir/out/err.txt",'a' )

with open(vcf_out + 'merged.vcf','a') as file_merged:
    bcftools_run = subprocess.Popen([bcftools,'merge','--merge','all','--force-samples',*vcf_path_lst_gz],universal_newlines=True,stdout=file_merged)#,stderr=logfile)
    bcftools_run.wait()


logfile.close()
print('END PIPLINE')