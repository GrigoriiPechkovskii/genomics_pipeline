import re
import sys
import csv
import os
import warnings

#import importlib

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt


print('start vcfproc')



class VcfData(pd.DataFrame):
    '''Main class for vcf (variant calling format) procsising '''

    def __init__(self,DataFrame=pd.DataFrame()):
        super().__init__()
        self.vcf = DataFrame

        #self.index = self['POS']
        #self.DataFrame = DataFrame


class VcfData():
    '''Main class for vcf (variant calling format) procsising '''
    COLUMNS_STANDARD = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
    COLUMNS_STANDARD_LENGTH = len(COLUMNS_STANDARD)

    def __init__(self,
        DataFrame=pd.DataFrame(columns=COLUMNS_STANDARD),
        set_uniq_index=True, drop_duplicate=True):
        """ Constructor for VcfData object
        VcfData required pd.DataFrame with vcf (variant calling format) format
        set_uniq_index : changing pd.DataFrame index on index = CHROM + POS + uniq number with '_' as delimiter
        """
        
        self.__vcf = self.__check_correctness_vcf(DataFrame)

        self.__vcf_bin = self.__vcf.iloc[:,9:].copy()#!!!must remove 

        self.vcf_drop_duplicate(drop_duplicate)
        self.vcf_reindexing(set_uniq_index)
        

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

    @property
    def vcf(self):
            return self.__vcf

    @vcf.setter
    def vcf(self,DataFrame):
        self.__vcf = self.__check_correctness_vcf(DataFrame)

        self.__vcf_bin = self.__vcf.iloc[:,9:].copy() #!!!must remove 


    @vcf.deleter
    def vcf(self):
        del self.__vcf

    def __check_correctness_vcf(self,DataFrame):
        """ Check type vcf, header etc...
        Type vcf must be pd.DataFrame
        Columns name from position 0 to 9 must be  ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
        If vcf correct return DataFrame of vcf
        """ 
        correct_pass = []
        if isinstance(DataFrame,pd.DataFrame):
            correct_pass.append(True)
        else:
            raise TypeError("Vcf must be pd.DataFrame")

        if all(DataFrame.columns[0:9] == self.COLUMNS_STANDARD):
            correct_pass.append(True)
        else:
            raise TypeError("Vcf columns is incorrect")

        index_na_exist = DataFrame[DataFrame['ALT'].isna()].index.append(DataFrame[DataFrame['REF'].isna()].index)
        if not index_na_exist.empty:
            DataFrame.drop(index_na_exist,inplace=True)
            warnings.warn("Warning vcf have NA value in REF or ALT columns. Deleted row = " + str(len(index_na_exist)),stacklevel=2)

        index_contains_non_standard_nucleotide = DataFrame[DataFrame['REF'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X').fillna(False)].index
        index_contains_non_standard_nucleotide = index_contains_non_standard_nucleotide.append(DataFrame[DataFrame['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X').fillna(False)].index)
        if not index_contains_non_standard_nucleotide.empty:
            DataFrame.drop(index_contains_non_standard_nucleotide,inplace=True)
            warnings.warn("Warning vcf contains non standard nucleotide NA value in REF or ALT columns. Deleted row = " + str(len(index_contains_non_standard_nucleotide)),stacklevel=2)

        columns_na_sum = sum(DataFrame.isna().any(axis=1))
        if columns_na_sum != 0:
            warnings.warn("Warning vcf contains NA in " + str(columns_na_sum) + " columns" ,stacklevel=2)

        if any(correct_pass) == True:
            return DataFrame
        else:
            raise TypeError("Vcf is incorrect")

    def vcf_reindexing(self,set_uniq_index=True):
        """ Reindexing vcf data auto
        """
        if set_uniq_index:
            num_str = np.array(range(self.__vcf.shape[0]),dtype=str)
            self.__vcf.index = self.__vcf['#CHROM'] + '_' + self.__vcf['POS'].astype(str) + '_'  + num_str
        else:
            self.__vcf.index = self.__vcf['#CHROM'] + '_' + self.__vcf['POS'].astype(str)

    def vcf_drop_duplicate(self,drop_duplicate=True):

        #if drop_duplicate:
        index_duplicated = self.__vcf[self.__vcf.duplicated(keep=False)].index
        if not index_duplicated.empty:
            if drop_duplicate:
                warnings.warn('Warning vcf have duplicate! Duplicate deleted, Deleted row = ' + str(len(index_duplicated)),stacklevel=2)
                self.__vcf.drop(index_duplicated, inplace=True)
            else:
                warnings.warn('Warning vcf have duplicate! Duplicate not deleted, Duplicate row = ' + str(len(index_duplicated)),stacklevel=2)






    def _samples_variation_template_slicer(self,samples_variation_dict):
        """ Function get samples_variation_dict and return pd.DataFrame sorted slice variation 
        samples_variation_dict  = {"Name variation":["CHROM", POS]}
        return pd.DataFrame sorted slice variation with Index from "Name variation"
        """
        #Here used loop for right sort
        template_samples_variation = pd.DataFrame()
        for key, val in samples_variation_dict.items():
            chrom_value = val[0]
            pos_value = int(val[1])

            template_samples_variation_val = self.__vcf[(self.__vcf['POS'] == pos_value) & (self.__vcf['#CHROM'] == chrom_value)]
            
            if not template_samples_variation_val.empty:
                template_samples_variation_val.index  = [key]
                template_samples_variation = template_samples_variation.append(template_samples_variation_val)
            else:
                raise ValueError("Variation " + key + " not found")

        self.template_samples_variation = template_samples_variation.iloc[:,9:]
        
        return self.template_samples_variation


    def to_genotype(self,template_genotype,samples_variation_dict):#!
        '''NA'''
        #samples_variants_for_genotype
        self.template_samples_variation = self._samples_variation_template_slicer(samples_variation_dict)

        if isinstance(template_genotype,pd.DataFrame):
            self.template_genotype = template_genotype
        else:
            raise TypeError("template_genotype must be pd.DataFrame")
        
        self.genotype = pd.Series()
        
        for template_var in self.template_genotype:
            for sample_var in self.template_samples_variation:
                df = pd.DataFrame({"Value_template_genotype":self.template_genotype[template_var], 
                                   "Value_sample_genotype":self.template_samples_variation[sample_var]})#!
                if ((df['Value_sample_genotype']).astype(str) == '.').any():
                    self.genotype[sample_var] = np.nan                    
                elif (df['Value_template_genotype'].astype(int) == df['Value_sample_genotype'].astype(int)).all():#! 222               
                    self.genotype[sample_var] = template_var
                    
        self.genotype = self.genotype[self.__vcf.columns[:self.COLUMNS_STANDARD_LENGTH]] #sort   



if __name__  == "__main__":


    samples_variation_dict = { 
                    'A.Br.001': ['NC_007530', 182106],
                    'A.Br.002': ['NC_007530', 947760],
                    'A.Br.003': ['NC_007530', 1493280],
                    'A.Br.004': ['NC_007530', 3600786],
                    'A.Br.006': ['NC_007530', 162509],
                    'A.Br.007': ['NC_007530', 266439],
                    'A.Br.008': ['NC_007530', 3947375],
                    'A.Br.009': ['NC_007530', 2589947],
                    'B.Br.001': ['NC_007530', 1458558],
                    'B.Br.002': ['NC_007530', 1056740],
                    'B.Br.003': ['NC_007530', 1494392],
                    'B.Br.004': ['NC_007530', 69952],
                    'A/B.Br.001': ['NC_007530', 3698013]}

    #vcf_inst = pd.read_csv('test_vcf.vcf',sep='\t',header=5)
    #vcf_inst = VcfData(vcf_inst.copy())
    #del vcf_reader

    df_can_ert = pd.read_table('can_snp_ert.csv', sep=',', engine='python',header=0)
    df_can_ert.index = df_can_ert['canSNP lineage/group']
    df_can_ert = df_can_ert.iloc[:,8:21].T   

    vcf_inst.to_genotype(df_can_ert,samples_variation_dict)
    vcf_inst.genotype 
