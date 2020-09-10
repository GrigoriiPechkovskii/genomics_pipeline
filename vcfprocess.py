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






if __name__  == "__main__":

    vcf_inst = pd.read_csv('test_vcf.vcf',sep='\t',header=5)
    vcf_inst = VcfData(vcf_inst.copy())