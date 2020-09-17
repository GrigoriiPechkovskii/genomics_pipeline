import os
import warnings
import time

import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

print('start vcfproc')

def timer_decor(func):
    def wrapper(*arg,**kwargs):
        start_time = time.time()
        loop = 200
        for _ in range(loop):
            result = func(*arg)
        
        result_time = round((time.time() - start_time),2)
        print("Function",func.__name__, "spend time in",loop , "loop = ", result_time, "second")
        return result
    return wrapper


class VcfData():
    '''Main class for vcf (variant calling format) procsising '''
    COLUMNS_STANDARD = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
    COLUMNS_STANDARD_LENGTH = len(COLUMNS_STANDARD)

    def __init__(self,
        DataFrame=pd.DataFrame(columns=COLUMNS_STANDARD),
        set_uniq_index=True, drop_duplicate=True, delete_variation_with_non_standard_nucleotide=True,
        drop_variation_with_na_in_ref_and_alt=True):
        """ Constructor for VcfData object
        VcfData required pd.DataFrame with vcf (variant calling format) format
        set_uniq_index : changing pd.DataFrame index on index = CHROM + POS + uniq number with '_' as delimiter
        """
        
        self.__vcf = self.__check_correctness_vcf(DataFrame,delete_variation_with_non_standard_nucleotide,
                                                  drop_variation_with_na_in_ref_and_alt)

        self.vcf_drop_duplicate(drop_duplicate)
        self.vcf_uniq_reindexing(set_uniq_index)


    @property
    def vcf(self):
            return self.__vcf

    @vcf.setter
    def vcf(self,DataFrame):
        self.__vcf = self.__check_correctness_vcf(DataFrame)

    @vcf.deleter
    def vcf(self):
        del self.__vcf

    @property
    def vcf_bin(self):
        """ vcf without standart calumns """
        return self.__vcf.iloc[:,self.COLUMNS_STANDARD_LENGTH:]

    def __check_correctness_vcf(self,DataFrame,delete_variation_with_non_standard_nucleotide=True,
                                drop_variation_with_na_in_ref_and_alt=True):
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
            raise TypeError("Vcf columns is incorrect", DataFrame.columns[0:9])

        index_na_exist = DataFrame[DataFrame['ALT'].isna()].index.append(DataFrame[DataFrame['REF'].isna()].index)
        if not index_na_exist.empty:
            if drop_variation_with_na_in_ref_and_alt:
                DataFrame.drop(index_na_exist,inplace=True)
                warnings.warn("Warning vcf have NA value in REF or ALT columns. Deleted row = " + str(len(index_na_exist)),stacklevel=2)
            else:
                warnings.warn("Warning vcf have NA value in REF or ALT columns. NOT Deleted row = " + str(len(index_na_exist)),stacklevel=2)
        #
        index_contains_non_standard_nucleotide = DataFrame[DataFrame['REF'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X').fillna(False)].index
        index_contains_non_standard_nucleotide = index_contains_non_standard_nucleotide.append(DataFrame[DataFrame['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X').fillna(False)].index)
        if not index_contains_non_standard_nucleotide.empty:
            if delete_variation_with_non_standard_nucleotide:
                DataFrame.drop(index_contains_non_standard_nucleotide,inplace=True)
                warnings.warn("Warning vcf contains non standard nucleotide NA value in REF or ALT columns. Deleted row = " + str(len(index_contains_non_standard_nucleotide)),stacklevel=2)
            else:
                warnings.warn("Warning vcf contains non standard nucleotide NA value in REF or ALT columns. Row not deleted.  Row  = " + str(len(index_contains_non_standard_nucleotide)),stacklevel=2)

        
        columns_na_sum = sum(DataFrame.iloc[:,self.COLUMNS_STANDARD_LENGTH:].isna().any(axis=1))
        if columns_na_sum != 0:
            warnings.warn("Warning vcf sample contains NA in " + str(columns_na_sum) + " columns" ,stacklevel=2)

        if all(correct_pass):
            return DataFrame
        else:
            raise TypeError("Vcf is incorrect")

    def vcf_uniq_reindexing(self,set_uniq_index=True):
        """ Reindexing vcf data auto
        """
        if set_uniq_index:
            num_str = np.array(range(self.__vcf.shape[0]),dtype=str)
            self.__vcf.index = self.__vcf['#CHROM'] + '_' + self.__vcf['POS'].astype(str) + '_'  + num_str
        else:
            self.__vcf.index = self.__vcf['#CHROM'] + '_' + self.__vcf['POS'].astype(str)

    def vcf_drop_duplicate(self,drop_duplicate=True):
        """ Delete row with duplicat CHROM, POS columns """
        index_duplicated = self.__vcf[self.__vcf[['#CHROM',"POS"]].duplicated(keep=False)].index
        if not index_duplicated.empty:
            if drop_duplicate:
                warnings.warn('Warning vcf have duplicate! Duplicate deleted, Deleted row = ' + str(len(index_duplicated)),stacklevel=2)
                self.__vcf.drop(index_duplicated, inplace=True)
            else:
                warnings.warn('Warning vcf have duplicate! Duplicate not deleted, Duplicate row = ' + str(len(index_duplicated)),stacklevel=2)

    def samples_variation_template_slicer(self,samples_variation_dict):
        """ Function get samples_variation_dict and return pd.DataFrame sorted slice variation,
        samples_variation_dict  = {"Name variation":["CHROM", POS]},
        return named_variation - dict {"index_vcf": "Name variation"}
        """
        #Here used loop for right sort
        named_variation = {}
        template_samples_variation = pd.DataFrame()
        for key, val in samples_variation_dict.items():
            chrom_value = val[0]
            pos_value = int(val[1])

            template_samples_variation_val = self.__vcf[(self.__vcf['POS'] == pos_value) & (self.__vcf['#CHROM'] == chrom_value)]
            
            if not template_samples_variation_val.empty:
                named_variation[template_samples_variation_val.index[0]] = key                 
            else:
                raise ValueError("Variation " + key + " not found")
        
        return named_variation


    def genotype_on_variation(self,template_for_genotype,samples_variation_dict):#!
        '''Determination genotype base on particular variation according to getting template,

        template_for_genotype  - pd.DataFrame consist in index name of variation, in columns name genotype,
        values is binary combination of int which generate genotype

        samples_variation_dict  - pd.DataFrame returned _samples_variation_template_slicer()
        return named_variation - dict {"index_vcf": "Name variation"}
        '''        
        named_variation = self.samples_variation_template_slicer(samples_variation_dict)
        self.template_samples_variation = self.vcf_bin.loc[named_variation].rename(index=named_variation)

        if isinstance(template_for_genotype,pd.DataFrame):
            self.template_for_genotype = template_for_genotype
        else:
            raise TypeError("template_for_genotype must be pd.DataFrame")
        
        self.genotype = pd.Series(dtype='object')
        
        for template_var in self.template_for_genotype:
            for sample_var in self.template_samples_variation:
                df = pd.DataFrame({"Value_template_genotype":self.template_for_genotype[template_var], 
                                   "Value_sample_genotype":self.template_samples_variation[sample_var]})#!

                if ((df['Value_sample_genotype']).astype(str) == '.').any():
                    self.genotype[sample_var] = np.nan                    
                elif (df['Value_template_genotype'].astype(int) == df['Value_sample_genotype'].astype(int)).all():#! 222               
                    self.genotype[sample_var] = template_var
                    
        self.genotype = self.genotype[self.__vcf.columns[self.COLUMNS_STANDARD_LENGTH:]] #sort

        return named_variation

    def determine_locus(self,locus_dir, reindex_variation_on_locus = False):
        """ """

        locus_df = pd.read_csv(locus_dir)
        determine_locus_index = []
        self.altname_variation = pd.Series(index = self.__vcf.index)
        for vcf_index in self.__vcf.index:
            position_variation = self.__vcf.loc[vcf_index]['POS']
            cotig_variation = self.__vcf.loc[vcf_index]['#CHROM']
            slice_locus = locus_df[(locus_df['start']<position_variation) & (position_variation<locus_df['end']) & (locus_df['contig']==cotig_variation)]
            if not slice_locus.empty:
                locus_find =  '_'.join(slice_locus['locus'])
                determine_locus_index.append(vcf_index)
                print(vcf_index, locus_find, vcf_index, position_variation, cotig_variation)
                self.altname_variation[vcf_index] = locus_find

        if reindex_variation_on_locus:
            determine_locus_index = self.vcf_altname_variation_reindexing()            

        return determine_locus_index
        
    def vcf_altname_variation_reindexing(self, delimiter="_"):
        """ """
        if hasattr(self, 'altname_variation'):            
            variation_with_altname =  self.altname_variation[self.altname_variation.notna()]            
            index_altname_dict = {index:index + delimiter + variation_with_altname[index] for index in variation_with_altname.index}
            self.__vcf.rename(index=index_altname_dict,inplace=True)

            return list(index_altname_dict.values())
        else:
            raise AttributeError("Vcf do not have altname_variation")

    def compute_param(self,number_variation=True, length_variation=True, mass_variation=True,delimiter="_"):
        """ """
        series_lst = []
        for num_row,variant_series in self.vcf.astype('object').iterrows():

            variant_lst = variant_series['REF'].split(',') + variant_series['ALT'].split(',')
            variant_dict = dict()
            for num, val in enumerate(variant_lst):
                param_lst = []
                if number_variation:
                    param_lst.append(str(num))
                if length_variation:
                    param_lst.append(str(len(val)))
                if mass_variation:
                    param_lst.append(str(self._lenmass(val)))

                variant_dict[str(num)] = delimiter.join(param_lst)
                variant_dict[num] = delimiter.join(param_lst) #!for int            

            variant_dict['.'] = '.'    
            series_lst.append(variant_series.astype('object').replace(variant_dict))

        self.vcf_param = pd.DataFrame(series_lst)


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


    def definition_core_vcf(self,type_core:str="TOTAL"):
        """  Definition core vcf from __vcf
            Core vcf - vcf do not have missing (NA or '.')values in sample columns
            Type core (str):  SNP - core have only SNP variation
                        TOTAL - core have all variation

        """
        self.core_vcf = self.__vcf[(self.vcf_bin.astype(str) != '.').all(axis=1)]
        if type_core == "TOTAL":            
            pass
        elif type_core == "SNP":
            self.core_vcf = self.core_vcf[self.core_vcf['INFO'].astype(str) == 'SNP']
        else:
            raise ValueError("Argument type_core must be 'TOTAL' or 'SNP'")

    def snp_uniq_finder(self):        
        
        dict_can_snp = samples_variation_dict

        if not hasattr(self, 'core_vcf'): 
            raise AttributeError("VcfData has no attribute 'core_vcf'")

        df_res = self.core_vcf.iloc[:,self.COLUMNS_STANDARD_LENGTH:].T
        
        snp_lst_uniq = []
        for i in df_res.columns:
            if i in df_res.columns:

                logic = df_res.eq(df_res[i], axis=0).all()#compare each snp with the entire dataframe
                a = logic[logic]

                if list(a.index.values) not in snp_lst_uniq:
                    snp_lst_uniq += [list(a.index.values)]
                    df_res.drop(a.index.values,axis=1,inplace=True)

        return snp_lst_uniq


    def to_set_snptype(self,named_set_snptype=None):        
        
        snp_lst_uniq = self.snp_uniq_finder()        
        self.snp_type = pd.Series(data = np.nan,index =self.vcf_bin.index,name='snp_type')
        uniq_number = 0               
        for n_uniq in range(len(snp_lst_uniq)):
            
            #This needs to change 
            if named_set_snptype != None:
                flag = True
                for named_snptype_key in named_set_snptype:
                    if named_snptype_key in snp_lst_uniq[n_uniq]:
                        self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])] = named_set_snptype[named_snptype_key]
                        flag = False
                if flag:
                    self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])] = self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])].replace(np.nan,'snp' + str(uniq_number+1))
                    uniq_number += 1
                    flag = True
            else:
                self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])] = self.snp_type[self.snp_type.index.isin(snp_lst_uniq[n_uniq])].replace(np.nan,'snp' + str(uniq_number+1))
                uniq_number += 1


def recluster_variant(vcf_data,variant_index,distance_pair=10):
    """ Functioin change value in vcf sample depending clustering (in present depending mass clustering)
        Input VcfData type,
        variant_index - vcf index variant which will be clustering
        distance_pair - min value for different cluster
        return clustered DataFrame 

    """
    
    new_cluster_df = pd.DataFrame(columns=vcf_data.vcf_bin.columns)
    
    for variant_index in variant_index:
        
        variation_condition_without_na = vcf_data.vcf.loc[variant_index][9:] != '.'
        slice_variant_len = vcf_data.vcf_param.loc[variant_index][9:][variation_condition_without_na]
        old_cluster_series = vcf_data.vcf.loc[variant_index][9:][variation_condition_without_na]

        slice_variant_len_reshape = np.reshape(slice_variant_len.values.astype(int), (len(slice_variant_len.values.astype(int)), 1))
        link = linkage(slice_variant_len_reshape,method='average', metric='euclidean')
        find_cluster = fcluster(link, distance_pair,criterion = 'distance')
        
        find_cluster_series = pd.Series(find_cluster,slice_variant_len.index)
        new_cluster_series = old_cluster_series.copy()

        zero_new = find_cluster_series.loc[old_cluster_series[old_cluster_series == '0'].index].values[0]#!check all values

        for i,j in zip(old_cluster_series.index,find_cluster_series.values):
            if new_cluster_series[i] == 0:
                continue
            elif j == zero_new:
                new_cluster_series[i] = 0
            else:
                new_cluster_series[i] = j

        new_cluster_series.update(new_cluster_series[new_cluster_series>zero_new] - 1)
        new_cluster_df = new_cluster_df.append(new_cluster_series)

    return new_cluster_df

def replace_non_standart_nucleotide_to_na_values(data):
        """ Function replace non standart nucleotide to na values (.)
            If variation cosnsist only  of variont with non standart nucleotide
            variation will delete
            Get data - pd.DataFrame like VcfData
            or VcfData(delete_variation_with_non_standard_nucleotide = False)
        """

        if isinstance(data,VcfData):
            data = data.vcf

        data_for_change = data[data['ALT'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X') |
                                      data['REF'].str.contains('N|W|R|M|Y|K|S|H|V|B|D|X')]
        data_for_change = data_for_change['REF'].str.split(',')  + data_for_change['ALT'].str.split(',')

        for var_lst,var_index in zip(data_for_change,data_for_change.index):
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
                    sclice_var_strain = data.loc[var_index][9:]
                    for var_strain, strain in zip(sclice_var_strain, sclice_var_strain.index):
                        if var_strain == '.':
                            continue
                        elif int(var_strain) > var_num:
                            data.loc[var_index, strain] = int(data.loc[var_index, strain]) - 1
                        elif int(var_strain) == var_num:
                            data.loc[var_index, strain] = '.'
                else:
                    new_alt.append(var_lst[var_num])                    
                data.loc[var_index, 'ALT'] = ','.join(new_alt[1:])
        data.drop(data[data['ALT'] == ''].index, inplace=True)
        data.drop(data[data['REF'] == ''].index, inplace=True)

        return data



if __name__  == "__main__":

    import unittest

    class TestVcfData(unittest.TestCase):

        def test_genotype(self):

            test_vcf_file = os.path.join('.','test','test_vcf.vcf')
            vcf_inst = pd.read_csv(test_vcf_file,sep='\t',header=5)
            vcf_inst = VcfData(vcf_inst.copy())

            test_inst_genotype = pd.Series({"B":"genotype1","C":"genotype2",
                                            "D":np.nan,"E":"genotype1",
                                            "K":"genotype3","Z":np.nan})
            samples_variation_dict = { 
                    'SNP01': ['A', 3],
                    'SNP02': ['A', 401],
                    'SNP03': ['A', 405],}

            template_for_genotype = pd.DataFrame(data=[[0,0,0],[1,1,0],[1,0,1]],columns=["genotype1","genotype2","genotype3"],
                                index=["SNP01","SNP02","SNP03"])

            named_variation = vcf_inst.genotype_on_variation(template_for_genotype,samples_variation_dict)

            self.assertEqual(list(vcf_inst.genotype.values), list(test_inst_genotype.values), "Incorrect genotype")
            self.assertEqual(list(vcf_inst.genotype.index), list(test_inst_genotype.index), "Incorrect genotype")

        def test_recluster_variant(self):

            test_vcf_file = os.path.join('.','test','test_vcf.vcf')
            vcf_inst = pd.read_csv(test_vcf_file,sep='\t',header=5)
            vcf_inst = VcfData(vcf_inst.copy())
            vcf_inst.compute_param(number_variation=False, length_variation=True, mass_variation=False,delimiter="_")

            new_cluster_df = recluster_variant(vcf_inst,["A_481_21","A_583_22"],distance_pair=10)
            recluster_df_test = pd.DataFrame(data=[[0, 0, 1, 1, np.nan, 1],[1, 0, 1, np.nan, 0, 0]],
                                             columns=['B', 'C', 'D', 'E', 'K', 'Z'],
                                             index=['A_481_21', 'A_583_22'])
            self.assertTrue((new_cluster_df.fillna(-1) == recluster_df_test.fillna(-1)).all(axis=None), "Incorrect genotype recluster variant")


    pd.options.display.max_columns = 20    

    samples_variation_dict = { 
                    'SNP01': ['A', 3],
                    'SNP02': ['A', 401],
                    'SNP03': ['A', 405],}

    template_for_genotype = pd.DataFrame(data=[[0,0,0],[1,1,0],[1,0,1]],columns=["genotype1","genotype2","genotype3"],
                                        index=["SNP01","SNP02","SNP03"])

    test_vcf_file = os.path.join('.','test','test_vcf.vcf')
    vcf_inst = pd.read_csv(test_vcf_file,sep='\t',header=0)
    #vcf_inst = vcf_inst.astype("object")
    vcf_inst = VcfData(vcf_inst.copy(),drop_duplicate=True,delete_variation_with_non_standard_nucleotide = True,
        drop_variation_with_na_in_ref_and_alt=True)
    
    '''
    template_for_genotype = pd.DataFrame(data=[[0,0,0],[1,1,0],[1,0,1]],columns=["genotype1","genotype2","genotype3"],
                                        index=["SNP01","SNP02","SNP03"])

    named_variation = vcf_inst.samples_variation_template_slicer(samples_variation_dict)
    named_variation = vcf_inst.genotype_on_variation(template_for_genotype,samples_variation_dict)
    vcf_inst.genotype

    vcf_inst.compute_param(number_variation=False, length_variation=True, mass_variation=False,delimiter="_")
    vcf_inst.vcf_param
    
 
    dv = replace_non_standart_nucleotide_to_na_values(vcf_inst)
    '''

    unittest.main(exit=False)