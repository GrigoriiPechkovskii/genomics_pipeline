import re
import sys
import csv
import os
import importlib

#import seaborn as sns
import numpy as np
import pandas as pd
#import matplotlib
import matplotlib.pyplot as plt
#from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
#import scipy

print('import vcf_process12')


class vcf_data():
    #def __init__(self,vcf):
    #   self.vcf = vcf
    #   self.vcf_bin = self.vcf.iloc[:,9:].copy()
    
    def __init__(self,vcf_dir,with_open=True,header_vcf=5,rename=True,dup_del=True,uniq_index=False):
        
        if with_open:
            self.vcf = pd.read_table(vcf_dir, sep='\t', engine='python',header=header_vcf)
        else:
            self.vcf = vcf_dir
        
        if uniq_index:
            num_str = np.array(range(self.vcf.shape[0]),dtype=str)
            self.vcf.index = self.vcf['#CHROM'] + '_' + self.vcf['POS'].astype(str) + '_'  + num_str
        else:
            self.vcf.index = self.vcf['#CHROM'] + '_' + self.vcf['POS'].astype(str)  
        
        
        self.vcf.drop(self.vcf[self.vcf['ALT'].isna()].index, inplace=True)#!!!!!
        
        #self.vcf = self.vcf[~self.vcf['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X')] #!!! check and del non standart nuc !!!!ERROR if locus have non standard nuc
        self.vcf = self.vcf[~self.vcf['REF'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X')] #!!! check and del non standart nuc
        



        data_for_change = self.vcf[self.vcf['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X')]
        data_for_change = data_for_change['REF'].str.split(',')  + data_for_change['ALT'].str.split(',')

        for var_lst,var_index in zip(data_for_change,data_for_change.index):
            #var_lst
            #{var_lst[i]:i for i in range(len(var_lst))}
            new_alt = []
            mem = 0
            for var_num in range(len(var_lst)):
                
                logic_non_standart = ("N" in var_lst[var_num] or "W" in var_lst[var_num] or
                                      "R" in var_lst[var_num] or "M" in var_lst[var_num] or
                                      "Y" in var_lst[var_num] or "K" in var_lst[var_num] or
                                      "S" in var_lst[var_num] or "H" in var_lst[var_num] or
                                      "V" in var_lst[var_num] or "B" in var_lst[var_num] or
                                      "D" in var_lst[var_num] or "X" in var_lst[var_num])        
                if logic_non_standart:

                    var_num -= mem
                    mem += 1                    
                    
                    sclice_var_strain = self.vcf.loc[var_index][9:]
                    
                    for var_strain, strain in zip(sclice_var_strain, sclice_var_strain.index):
                        if var_strain == '.':
                            continue
                        elif int(var_strain) > var_num:
                            self.vcf.loc[var_index, strain] = int(self.vcf.loc[var_index, strain]) - 1
                        elif int(var_strain) == var_num:
                            self.vcf.loc[var_index, strain] = '.'
                else:
                    new_alt.append(var_lst[var_num])            
                    
                    
                self.vcf.loc[var_index, 'ALT'] = ','.join(new_alt[1:])

        self.vcf.drop(self.vcf[self.vcf['ALT'] == ''].index, inplace=True)




        #self.vcf.index = self.vcf['POS']
        self.vcf.columns = [i.replace('.fna','').replace('.gbk','') for i in list(self.vcf.columns)]
        #check NA
        print('Checking NA =',sum(self.vcf.isna().any()))
        
        if rename == True:
            #def rename_assemly(self):        
            table_read = pd.read_csv('genomes_proks.csv')
            #assemb = pd.DataFrame(list(table_read['Assembly'].str.strip().str.replace('GCA','GCF').str.split('.').values))[0].values
            assemb = pd.DataFrame(list(table_read['Assembly'].str.strip().str.replace('GCA','GCF').values))
            table_read['Strain'] = table_read['Strain'].str.strip().str.replace(' ','_') + '_' + assemb[0]

            rep = table_read[['Strain']]
            rep.index = assemb[0]
            strain_dict = rep['Strain'].to_dict()

            q = pd.DataFrame(data = list( self.vcf.columns.str.split('_',2).values))
            q = q[0] +'_'+ q[1]
            q.dropna(inplace=True)
            #w = pd.Series(strain_dict)[q].values
            e = list(self.vcf.columns[:9]) + list(pd.Series(strain_dict)[q].values)
            self.vcf.columns = e
        
        
        
        self.vcf_bin = self.vcf.iloc[:,9:].copy()
        
        #df_vcf.index.get_duplicates()
        if dup_del:
            ind_dup = self.vcf.index[self.vcf.index.duplicated()].unique()
            if len(ind_dup) != 0:
                print('Warning vcf have duplicate! Duplicate deleted',len(ind_dup))
                self.vcf = self.vcf.drop(ind_dup)            
            
        self.dict_can_snp = {'NC_007530_182106':'A.Br.001',
                'NC_007530_947760':'A.Br.002',
                'NC_007530_1493280':'A.Br.003',
                'NC_007530_3600786':'A.Br.004',
                'NC_007530_162509':'A.Br.006',
                'NC_007530_266439':'A.Br.007',
                'NC_007530_3947375':'A.Br.008',
                'NC_007530_2589947':'A.Br.009',
                'NC_007530_1458558':'B.Br.001',#!
                'NC_007530_1056740':'B.Br.002',
                'NC_007530_1494392':'B.Br.003',
                'NC_007530_69952':'B.Br.004',
                'NC_007530_3698013':'A/B.Br.001'}

        if False:
            self.dict_can_snp = {'NC_007530_182106':'A.Br.001',
                    'NC_007530_947760':'A.Br.002',
                    'NC_007530_1493280':'A.Br.003',
                    'NC_007530_3600786':'A.Br.004',
                    'NC_007530_162509':'A.Br.006',
                    'NC_007530_266439':'A.Br.007',
                    'NC_007530_3947375':'A.Br.008',
                    'NC_007530_2589947':'A.Br.009',
                    'NC_007530_1458558':'B.Br.001',#!
                    'NC_007530_1056740':'B.Br.002',
                    'NC_007530_1494392':'B.Br.003',
                    'NC_007530_69952':'B.Br.004',
                    'NC_007530_3698013':'A/B.Br.001'}

        self.dict_can_snp2 = {'NC_007530_182106_806_A.Br.001':'A.Br.001',
                'NC_007530_947760_4824_A.Br.002':'A.Br.002',
                'NC_007530_1493229_7411_A.Br.003':'A.Br.003',
                'NC_007530_3600767_22240_A.Br.004':'A.Br.004',
                'NC_007530_162509_712_A.Br.006':'A.Br.006',
                'NC_007530_266439_1238_A.Br.007':'A.Br.007',
                'NC_007530_3947375_23703_A.Br.008':'A.Br.008',
                'NC_007530_2589947_15116_A.Br.009':'A.Br.009',
                'NC_007530_1455386_7265_B.Br.001':'B.Br.001',#!
                'NC_007530_1056740_5354_B.Br.002':'B.Br.002',
                'NC_007530_1494392_7415_B.Br.003':'B.Br.003',
                'NC_007530_69952_344_B.Br.004':'B.Br.004',
                'NC_007530_3697997_22562_A/B.Br.001':'A/B.Br.001'}

        
    def rename_assemly(self):
        
        table_read = pd.read_csv('genomes_proks.csv')
        #assemb = pd.DataFrame(list(table_read['Assembly'].str.strip().str.replace('GCA','GCF').str.split('.').values))[0].values
        assemb = pd.DataFrame(list(table_read['Assembly'].str.strip().str.replace('GCA','GCF').values))
        table_read['Strain'] = table_read['Strain'].str.strip().str.replace(' ','_') + '_' + assemb[0]

        rep = table_read[['Strain']]
        rep.index = assemb[0]
        strain_dict = rep['Strain'].to_dict()

        q = pd.DataFrame(data = list( self.vcf.columns.str.split('_',2).values))
        q = q[0] +'_'+ q[1]
        q.dropna(inplace=True)
        #w = pd.Series(strain_dict)[q].values
        e = list(self.vcf.columns[:9]) + list(pd.Series(strain_dict)[q].values)
        self.vcf.columns = e

    
    def can_slice(self):        
        #Sorting automatically
        #Slice for canonical SNPs
        
        #df_slice_can2 =  self.vcf.loc[['NC_007530_182106','NC_007530_947760','NC_007530_1493280',
        #                  'NC_007530_3600786','NC_007530_162509','NC_007530_266439',
        #                  'NC_007530_3947375','NC_007530_2589947','NC_007530_1458558',
        #                 'NC_007530_1056740','NC_007530_1494392','NC_007530_69952','NC_007530_3698013']]
        df_slice_can = pd.DataFrame()
        for key in self.dict_can_snp.keys():
            key_split = key.split('_')
            chrom_key = key_split[0]+'_'+key_split[1]
            pos_key = int(key_split[2])

            df_slice_can = df_slice_can.append(self.vcf[(self.vcf['POS'] == pos_key) & (self.vcf['#CHROM'] == chrom_key)])

        #df_slice_can =  self.vcf.loc[list(self.dict_can_snp.keys())]
        #print(df_slice_can)
        self.df_can = df_slice_can.iloc[:,9:]
        
        #if True:
        #    self.df_can = self.df_can.astype(str).replace('2',1).astype(int)#! warning for if can sclice have severral variance

        #self.df_can = self.df_can.rename(index = self.dict_can_snp)
        self.df_can.index = list(self.dict_can_snp.values())
        
        return self.df_can

    #The presence of a unique snp in the core genome
    #df_res = df_vcf_extend.vcf.T.copy()
    def snp_uniq_finder(self,df_res):
        
        if True:
            dict_can_snp = self.dict_can_snp
        
        snp_lst_uniq = []
        for i in df_res.columns:
            if i in df_res.columns:

                logic = df_res.eq(df_res[i], axis=0).all()#compare each snp with the entire dataframe
                a = logic[logic]#clice

                if list(a.index.values) not in snp_lst_uniq:
                    snp_lst_uniq += [list(a.index.values)]
                    df_res.drop(a.index.values,axis=1,inplace=True)

        #Determining the frequency of snp in the cow genome (in conjunction with canSnp)
        #The name of the group of identical SNPs was taken from the first position in this group
        freq_snp = {}
        snp_used = []
        count = 0
        for snp_uniq in snp_lst_uniq:

            for can_snp in dict_can_snp:#separately for canonically
                if can_snp in snp_uniq:
                    snp_used += [snp_uniq]
                    count += 1
                    n = len(snp_uniq)
                    freq_snp[dict_can_snp[can_snp]] = n
            if snp_uniq not in snp_used:
                count += 1
                n = len(snp_uniq)
                freq_snp[snp_uniq[0]] = n
        return snp_lst_uniq

            
    def genotype_on_locus(self,locus_dir): 
    
        locus_df = pd.read_csv(locus_dir)#!
        genotype_on_locus_index = []
        for vcf_index in self.vcf.index:
            locus_pos = self.vcf.loc[vcf_index]['POS']
            cotig = self.vcf.loc[vcf_index]['#CHROM']
            slice_vcf = locus_df[(locus_df['start']<locus_pos) & (locus_pos<locus_df['end']) & (locus_df['contig']==cotig)]
            if not slice_vcf.empty:
                locus_find =  '_'.join(slice_vcf['locus'])
                print(vcf_index + '_' + locus_find,vcf_index,locus_pos,cotig)
                self.vcf.rename(index={vcf_index:vcf_index + '_' + locus_find},inplace=True)            
                genotype_on_locus_index.append(vcf_index + '_' + locus_find) 
        return genotype_on_locus_index
    
    def to_genotype(self,template_genotype,sample_genotype=None):#!
        '''NA'''
        
        self.sample_genotype = sample_genotype
        if True:
            self.can_slice()
            self.sample_genotype = self.df_can

        self.template_genotype = template_genotype
        self.genotype = pd.Series()
        
        
        for template_var in self.template_genotype:
            for sample_var in self.sample_genotype:
                df = pd.DataFrame({"Value_template_genotype":self.template_genotype[template_var], 
                                   "Value_sample_genotype":self.sample_genotype[sample_var]})#!
                if ((df['Value_sample_genotype']).astype(str) == '.').any():
                    #print(df['Value_sample_genotype'])#222
                    continue


                if (df['Value_template_genotype'].astype(int) == df['Value_sample_genotype'].astype(int)).all():#! 222               
                    self.genotype[sample_var] = template_var
                    
        self.genotype = self.genotype[self.vcf_bin.columns] #sort   
    
    
    
    def to_set_snptype(self,dict_can_snp=None):
        
        if True:
            dict_can_snp = self.dict_can_snp
        
        snp_lst_uniq = self.snp_uniq_finder(self.vcf_bin.T)
        
        self.snp_type = pd.Series(data = 'unknown',index =self.vcf_bin.index,name='snp_type')        
               
        for n_uniq in range(len(snp_lst_uniq)):
            #condition for a canonical label
            if len(set(dict_can_snp.keys()) & set(snp_lst_uniq[n_uniq])) > 0:
                label = dict_can_snp[list(dict_can_snp.keys() & set(snp_lst_uniq[n_uniq]))[0]]
                self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])] = self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])].replace('unknown',label)

            else:
                self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])] = self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])].replace('unknown','snp' + str(n_uniq+1))
                #self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])]
                
    def be_snp_genome_counter(self):
        
        self.be_snp_genome_count = pd.Series(self.vcf_bin.astype(float).sum(axis=1),name='be_snp_genome_count' )        
        
    def _lenvar(self,var):
        if pd.notna(var):
            return len(var)
        else:
            return None

    def _lenmass(self,var):
        nuc_mass = {'A':331.2,'C':307.2,'G':347.2,'T':322.2}
        mass = 0    
        if pd.notna(var):
            for v in var:
                mass += nuc_mass[str(v)]
            return round(mass,2)
        else:
            return None

    def compute_lenmass(self):
        df = pd.DataFrame(data = list(self.vcf['ALT'].str.split(',').values),index = self.vcf.index,)
        df.fillna(value=pd.np.nan, inplace=True)
        alt_name = []
        for n in range(1,df.shape[1]+1):
                  alt_name += ['alt_seq_' + str(n)]                               
        df.columns = alt_name
        df['alt_seq_ref_0'] = self.vcf['REF'].values
        cols = df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df = df[cols]

        df_len = pd.DataFrame()
        df_mass = pd.DataFrame()
        for num,col in enumerate(df):
            mass_name = 'alt_mass_' + str(num)
            len_name = 'alt_len_' + str(num)
            df_len[len_name] = df[col].apply(self._lenvar)
            df_mass[mass_name] = df[col].apply(self._lenmass)
        

        df_len_max = pd.DataFrame([(df_len['alt_len_0'] - df_len[i]).abs() for i in df_len.iloc[:,1:]]).T.max(axis=1)
        df_len_max.name = 'alt_len_diff_max'

        df_mass_max = pd.DataFrame([(df_mass['alt_mass_0'] - df_mass[i]).abs() for i in df_mass.iloc[:,1:]]).T.max(axis=1)
        df_mass_max.name = 'alt_mass_diff_max'

        self.vcf_lenmass = pd.concat([df,df_len,df_len_max,df_mass,df_mass_max],axis=1,sort=False)
        self.df_len = df_len
        self.df_mass = df_mass
        
        self.seq_variance = df
        #return self.vcf_lenmass
        
#df_vcf_extend = df_vcf_extend.rename(index = dict_can_snp)

    def compute_binlen (self):
        #self.vcf_binlen = pd.DataFrame()
        series_lst = []
        for num,variant_series in self.vcf.astype('str').iterrows():
            variant_lst = variant_series['REF'].split(',') + variant_series['ALT'].split(',')
            variant_dict = {str(num):str(num) +'_'+ str(len(val)) + '_' + str(self._lenmass(val)) for num,val in enumerate(variant_lst)}
            variant_dict['.'] = '.'    
            series_lst.append(variant_series.replace(variant_dict))
        self.vcf_binlen = pd.DataFrame(series_lst)
        
        self.vcf_varlen = pd.DataFrame(series_lst)
        
    def compute_varlen (self):
        #self.vcf_binlen = pd.DataFrame()
        series_lst = []
        for num,variant_series in self.vcf.astype(str).iterrows():
            variant_lst = variant_series['REF'].split(',') + variant_series['ALT'].split(',')
            variant_dict = {str(num):len(val) for num,val in enumerate(variant_lst)}
            variant_dict['.'] = '.'    
            series_lst.append(variant_series.astype(str).replace(variant_dict))

        self.vcf_varlen = pd.DataFrame(series_lst)
        
        
    def variance_annotation_processing (self):
        ''' '''
        ann_name = 'Allele|Annotation|Annotation_Impact|Gene_Name|id|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos / cDNA.length|CDS.pos / CDS.length|AA.pos / AA.length|Distance|ERRORS / WARNINGS / INFO'.split('|') 
        #print('1111')
        type_ann_type = pd.DataFrame(data = list(self.vcf['INFO'].str.split(';ANN=',1).values),index = self.vcf.index)
        #print(type_ann_type)
        type_ann_type[0].name = 'type_variance'
        self.type_variance = type_ann_type[0]
        type_ann_type[1].fillna('inv',inplace=True)
        self.type_ann = pd.DataFrame(data = list(type_ann_type[1].str.split(',').values),index = self.vcf.index)
        #self.type_ann = type_ann
        self.type_ann[0] = self.type_ann[0].replace({'inv':np.nan})

        self.annotation_lst = []
        self.name_annotation = pd.DataFrame(index=self.vcf.index)
        self.id_annotation = pd.DataFrame(index=self.vcf.index)
        self.allele_annotation = pd.DataFrame(index=self.vcf.index)
        for i in range(self.type_ann.shape[1]):    
            ann_slice = pd.DataFrame(data = list(self.type_ann[i].dropna().str.split('|').values),index = self.type_ann[i].dropna().index).iloc[:,:16]
            if ann_slice.shape[1] != len(ann_name):
                print('Warning not eq len in annotation, will be pass this annotation!',ann_slice)
                continue
            else:
                ann_slice.columns = ann_name
                self.name_annotation['Annotation_' + str(i)] = ann_slice['Annotation']
                self.id_annotation ['id_' + str(i)] = ann_slice['id']
                self.allele_annotation ['Allele_' + str(i+1)] = ann_slice['Allele'].str.replace('ANN=','')
                self.annotation_lst += [ann_slice]

        '''self.seq_variance
        print(all(self.name_annotation.index == self.seq_variance.index))
        for index in self.name_annotation.index:
            if len(self.name_annotation.loc[index].dropna()) == len(self.seq_variance.loc[index].dropna()) - 1:
                pass'''
            #else:
                #print('Warning count annotations and seq_variance not equal',index)

    
    def gff_merge(self,gff_dir):

        df_annot_gff = pd.read_table(gff_dir, engine='python',skiprows=9, skipfooter=1,names=['contig','source','type_seq','start','end','5','strand','7','features'])
        #удаление строк с #
        df_annot_gff.drop(df_annot_gff[df_annot_gff['contig'].str.contains('#')].index,inplace=True)
        index_start_end = df_annot_gff['start'].astype('int').astype('str') + '_to_' + df_annot_gff['end'].astype('int').astype('str')
        df_annot_gff.index = index_start_end


        #!!!LONG time
        #the creation of a dictionary - type of sequence : type of feature
        df_type_seq = pd.DataFrame()

        type_seq_features = dict()
        type_seq_set = set(df_annot_gff['type_seq'])

        df_seq_var = pd.DataFrame()

        #checking features just in case for each type of sequence
        name_features = set()
        for line,ind in zip(df_annot_gff['features'],df_annot_gff['features'].index):
            result = ['id'] + re.findall(r';(.*?)=', line)
            result = [i.lower() for i in result]


            s = pd.Series(data=line.split(';'),index=result,name=ind)
            s = s[~s.index.isin(['is_circular','genome','mol_type','strain','old_locus_tag',
                                 'parent','inference','transl_table','anticodon',
                                'end_range','partial','start_range','bound_moiety','regulatory_class',
                                'gene_synonym','exception','ribosomal_slippage','old-name','plasmid-name'])]
            df_seq_var = df_seq_var.append([s],sort=False)
        #combining a gff piece to features and a piece with separated features
        df_annot_gff_full =pd.concat([df_annot_gff.loc[:,:'7'],df_seq_var],axis=1, sort=False)
        gff_id = pd.DataFrame(data=list(df_annot_gff_full['id'].str.split('-').values),index=df_annot_gff_full.index)[1]#replace('ID=',''))
        df_annot_gff_full['id'] = gff_id


        self.annotation_lst[0]
        result = pd.merge(self.annotation_lst[0], df_annot_gff_full[(df_annot_gff_full['type_seq'] == 'gene') | (df_annot_gff_full['type_seq'] == 'pseudogene')],
                          how='left',on='id',left_index=False)

        result.index = self.annotation_lst[0].index#only 1 list !!!!!!!!!!!!!!!!!
        
        self.ann_merge = result
        #return result
        
    def standart_processing(self):
        
        #self.to_genotype(df_can_ert)#,df_can)
        self.to_set_snptype()#dict_can_snp)
        #self.be_snp_genome_counter()
        self.compute_lenmass()
        #self.variance_annotation_processing()

        gff_dir = 'test_GCF_000008445.1_ASM844v1_genomic.gff'
        #self.gff_merge(gff_dir)

        
def recluster_variant(vcf_data,variant_index,distance_pair=10):
    
    new_cluster_df = pd.DataFrame()
    new_cluster_df = pd.DataFrame(columns=vcf_data.vcf_bin.columns)

    
    for variant_index in variant_index:
    
        slice_variant_len = vcf_data.vcf_varlen.loc[variant_index][9:][vcf_data.vcf_varlen.loc[variant_index][9:] != '.']

        slice_variant_bin_old = vcf_data.vcf.loc[variant_index][9:][vcf_data.vcf.loc[variant_index][9:] != '.']
        slice_variant_undet = vcf_data.vcf_varlen.loc[variant_index][9:][vcf_data.vcf_varlen.loc[variant_index][9:] == '.']
        #slice_variant_zero = vcf_data.vcf_bin.loc[variant_index][9:][vcf_data.vcf_bin.loc[variant_index][9:] == '0']

        slice_variant_len_reshape = np.reshape(slice_variant_len.values.astype(int), (len(slice_variant_len.values.astype(int)), 1))

        link = linkage(slice_variant_len_reshape,method='average', metric='euclidean')
        find_cluster = fcluster(link, distance_pair,criterion = 'distance')
        
        find_cluster_series = pd.Series(find_cluster,slice_variant_len.index)
        old_cluster_series = slice_variant_bin_old        
        
        new_cluster_series = slice_variant_bin_old.copy()        
        
        #zero_new = find_cluster_series.iloc[old_cluster_series[old_cluster_series == 0]].values[0]#!check all values
        zero_new = find_cluster_series.loc[old_cluster_series[old_cluster_series == 0].index].values[0]#!check all values

        for i,j in zip(old_cluster_series.index,find_cluster_series.values):
            if new_cluster_series[i] == 0: #or slice_NC_007530_226253_1006_bin[i] == '.':
                #print(i)
                continue
            elif j == zero_new:
                new_cluster_series[i] = 0
            else:
                new_cluster_series[i] = j           
        
        new_cluster_series.update(new_cluster_series[new_cluster_series>zero_new] - 1)
        
       
        new_cluster_df = new_cluster_df.append(new_cluster_series)       
        
        
    
    return new_cluster_df




def recluster_variant(vcf_data,variant_index,distance_pair=10):
    
    new_cluster_df = pd.DataFrame()
    new_cluster_df = pd.DataFrame(columns=vcf_data.vcf_bin.columns)

    
    for variant_index in variant_index:
    
        slice_variant_len = vcf_data.vcf_varlen.loc[variant_index][9:][vcf_data.vcf_varlen.loc[variant_index][9:] != '.']

        slice_variant_bin_old = vcf_data.vcf.loc[variant_index][9:][vcf_data.vcf.loc[variant_index][9:] != '.']
        slice_variant_undet = vcf_data.vcf_varlen.loc[variant_index][9:][vcf_data.vcf_varlen.loc[variant_index][9:] == '.']
        #slice_variant_zero = vcf_data.vcf_bin.loc[variant_index][9:][vcf_data.vcf_bin.loc[variant_index][9:] == '0']

        slice_variant_len_reshape = np.reshape(slice_variant_len.values.astype(int), (len(slice_variant_len.values.astype(int)), 1))

        link = linkage(slice_variant_len_reshape,method='average', metric='euclidean')
        find_cluster = fcluster(link, distance_pair,criterion = 'distance')
        
        find_cluster_series = pd.Series(find_cluster,slice_variant_len.index)
        old_cluster_series = slice_variant_bin_old        
        
        new_cluster_series = slice_variant_bin_old.copy()        
        
        #zero_new = find_cluster_series.iloc[old_cluster_series[old_cluster_series == 0]].values[0]#!check all values
        zero_new = find_cluster_series.loc[old_cluster_series[old_cluster_series == 0].index].values[0]#!check all values

        for i,j in zip(old_cluster_series.index,find_cluster_series.values):
            if new_cluster_series[i] == 0: #or slice_NC_007530_226253_1006_bin[i] == '.':
                #print(i)
                continue
            elif j == zero_new:
                new_cluster_series[i] = 0
            else:
                new_cluster_series[i] = j           
        
        new_cluster_series.update(new_cluster_series[new_cluster_series>zero_new] - 1)
        
       
        new_cluster_df = new_cluster_df.append(new_cluster_series)       
        
        
    
    return new_cluster_df
  
    
