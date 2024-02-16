#########################################################################
#                          IMPORTED MODULES                             #
#########################################################################

from os import lstat
import pandas as pd
import numpy as np
import copy
import pickle
#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import StrVector
#from rpy2.robjects import pandas2ri
#import rpy2.robjects as robjects


#The following functions were originally written for the purposes of the:
#########################################################################
#                        naming_files2CCLE                              #
#########################################################################
#python script                                                    

def LOAD2DFCOLUMNS(file_path, columns):
    '''
    This function loads a dataframe and subsets specified columns.
    This function requires pandas module imported as follows:
    > import pandas as pd
    ----------------------------------------------------------------
    INPUT
         <<file_path>>: Path to the dataframe file (tab seperated file)
         <<columns>>: List of columns(s) to subset 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> interest: DataFrame containing only the wanted 
                               columns.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    df = pd.read_csv(file_path, sep = '\t')
    interest = df[columns]
    return interest

def SELECT_ROWS_BY_VALUE_EQUALITY(df, column_to_check, target):
    '''
    This function subsets rows based on the content of specific row
    in a specific column.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<df>>: Dataframe
         <<column_to_check>>: Name of the column to iterate over.
         <<target>>: Value/String used to select rows. 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> interest: DataFrame containing only the wanted 
                               rows.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    interest = df[df[column_to_check] == target]
    return interest

def INDEX_OF_ROWS_CONTAINING_SUBSTRING(df, column_to_check, target):
    '''
    This function gives indexes of rows containing a specific substring
    in a specific column.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<df>>: Dataframe
         <<column_to_check>>: Name of the column to iterate over.
         <<target>>: Value/String used to select rows. 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> idx: List of indexes of rows containing substring.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    res = df.index[df[column_to_check].str.contains(target)]
    cpt = len(res)
    idx = []
    if cpt != 0: 
        while cpt >= 0:
            idx.append(res[cpt-1])
            cpt -= 1
    return idx 

def DROP_ROWS_BY_INDEX(df, idx):
    '''
    This function deletes rows using dataframe indexes.
    This function requires pandas module imported as follows:
    > import pandas as pd
    ----------------------------------------------------------------
    INPUT
         <<df>>: DataFrame
         <<columns>>: List of indexes 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> shaved: DataFrame containing only the wanted 
                             rows.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    
    if len(idx) != 0:
        shaved = df.drop(idx)
    else:
        shaved = df
    return shaved

def TRIM_FILE_NAME(chain, front, back):
    '''
    This function trims the extremities of a string.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<chain>>: String to trim
         <<front>>: Index of the last character trimmed in the front.
         <<back>>: Index of the last character trimmed in the back 
                   (if negative index is used). 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> new: Trimmed string.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    new = chain[front : back]
    return new


#The following functions were originally written for the purposes of the:
#########################################################################
#                         sampleContext2R.py                            #
#########################################################################

def FILEREADER2LIST(path_list):
    '''
    This function reads text file(s) (one element per line) into a list.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<path_list>>: List of strings containing path(s) to the text
                        files
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> lst: List of elements in all the text files.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    september 2, 2022
    '''
    lst = []
    for path in path_list:
        with open(path) as file:
            for line in file:
                lst.append(line.strip())   
    return lst

def COLUMNSELECTION(df, lst):
    '''
    This function selects specific columns of a dataframe.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<df>>: Dataframe to be subset.
         <<lst>>: List of columns used to subset
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         sub -> DataFrame: df subset
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    september 2, 2022
    '''
    sub = df[lst]
    return sub

#The following functions were originally written for the purposes of the:
#########################################################################
#                         contextOrdering.py                            #
#########################################################################

def CONTEXTEORDERING(order, mut_mat):
    '''
    This function reorders rows of a dataframe.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
	<<order>>: File path to row order.
	<<mut_mat>>: Dataframe to reorder
	________________________________________________________________

	-----------------------------------------------------------------
	RETURNS
	ordered -> DataFrame: reordered dataframe.
	_________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    september 2, 2022
	'''
    lst = [order]
    sequence = FILEREADER2LIST(lst)
    df = pd.read_csv(mut_mat, sep='\t')
    idf = df.set_index('MutationType')
    ordered = idf.reindex(sequence)
    ordered.insert(loc=0, column='MutationType', value=sequence)
    return ordered

#The following functions were originally written for the purposes of the:
#########################################################################
#                         XXXXXXXXXXXXXXX.py                            #
#########################################################################

def TISSUE_PARSER(dependency, tissue):
    '''
    This function selects columns of samples of a specific name.
    (based on CCLE names)
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
	<<dependency>>: DataFrame: Column names are sample CCLE names.
	<<tissue>>: String of the tissue type to select (tissue type 
                determined by nik-zainal).
	________________________________________________________________

	-----------------------------------------------------------------
	RETURNS
	df -> DataFrame:
	_________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    september 24, 2022
	'''
    tissues = {"Biliary":['BILIARY_TRACT'], "Bladder":['URINARY_TRACT'], "Bone_SoftTissue":['SOFT_TISSUE', 'BONE'], "Breast":['BREAST'], "CNS":['AUTONOMIC_GANGLIA', 'CENTRAL_NERVOUS_SYSTEM'], "Colorectal":['LARGE_INTESTINE'], "Esophagus":['OESOPHAGUS'], "Head_neck":[], "Kidney":['KIDNEY'], "Liver":['LIVER'], "Lung":['LUNG', 'PLEURA'], "Lymphoid":['HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'], "Myeloid":['HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'], "NET":['THYROID'], "Oral_Oropharyngeal":['SALIVARY_GLAND', 'UPPER_AERODIGESTIVE_TRACT'], "Ovary":['OVARY'], "Pancreas":['PANCREAS'], "Prostate":['PROSTATE'], "Skin":['SKIN'], "Stomach":['STOMACH'], "Uterus":['CERVIX', 'ENDOMETRIUM'], "Other":['PLACENTA']}
    cols = ['MutationType']
    t = tissues[tissue]
    for i in t:
        for c in dependency:
            if i in c:
                cols.append(c)
    df = dependency[cols]
    return df

#The following functions were originally written for the purposes of the:
#########################################################################
#                       object2vcf(cluster).py                          #
#########################################################################
def MERGE_GROUPED_MUTATIONS(df):
    df = df.reset_index(drop=True)
    lines = df.index.to_list()
    limit = lines[-1]
    refind = df.columns.to_list().index('REF')
    altind = df.columns.to_list().index('ALT')
    posind = df.columns.to_list().index('POS')
    to_skip = []
    keep = []
    for ind in lines:
        if ind not in to_skip:
            ref = df.iloc[ind, refind]
            alt = df.iloc[ind, altind]
            pos = df.iloc[ind, posind]
            if ind+1 <= limit and int(pos)+1 == int(df.iloc[ind+1, posind]):
                keep.append(ind)
                n = 1
                while ind+n <= limit and int(pos)+n == int(df.iloc[ind+n, posind]):
                    to_skip.append(ind+n)
                    ref = ref + str(df.iloc[ind+n, refind])
                    alt = alt + str(df.iloc[ind+n, altind])
                    n += 1
            else:
                keep.append(ind)
            df.iloc[ind, refind] = ref
            df.iloc[ind, altind] = alt
    gbs = df.loc[keep]
    return gbs        


def DATA_CLEANER(vcf, output, sample):
    '''
    This function converts a vcf file to a dataframe while only
    keeping 4 columns.
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<vcf>>: vcf file path.
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	chr_pos_ref_alt_info -> DataFrame: Conserved columns: #CHROM, POS, REF, ALT.
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 7, 2022
	'''
    vcf_file =  open(vcf, 'r')
    data_file =  open(f'{vcf}.data', 'w')
    for line in vcf_file:
        if line.startswith('#CHROM'):
            data_file.write(line)
        if not (line.startswith('#')):
            data_file.write(line)
    vcf_file.close()
    data_file.close()
    all_mut_info = pd.read_csv(f'{vcf}.data', sep='\t')
    chr_pos_ref_alt_info = all_mut_info[['#CHROM', 'POS', 'REF', 'ALT']]
    chromo = set(chr_pos_ref_alt_info['#CHROM'])
    for chr in chromo:
        single_chrom = chr_pos_ref_alt_info[chr_pos_ref_alt_info['#CHROM']== chr]
        clean_single_chrom = MERGE_GROUPED_MUTATIONS(single_chrom)
        clean_single_chrom.to_csv(f'{output}/{sample}.data.chr{chr}', sep = '\t', index = False)
    return f'{vcf}.data'

def MID_GENERATOR(data):
    '''
    This function generates identifiers for all mutations contained 
    in a dataframe. These IDs were meant for the generation of objects
    produced by mutationClassification.py script.
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<data>>: Mutation dataframe produced from a vcf file 
    (see DATA_CLEANER function output).
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	MID_lst -> List of all IDs.
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 7, 2022
	'''
    vcf_data = pd.read_csv(data, sep = '\t')
    row = 0
    MID_lst = []
    for i in vcf_data['#CHROM']:
        MID_lst.append('MID.chr'+str(i)+'_'+str(vcf_data.loc[row][1])+'_'+vcf_data.loc[row][2]+'_'+vcf_data.loc[row][3])
        row += 1 
    return MID_lst 

def MID_GENERATOR_FROM_DATAFRAME(df):
    df = df.reset_index()
    row = 0
    MID_lst = []
    for i in df['#CHROM']:
        MID_lst.append('MID.chr'+str(i)+'_'+str(df.loc[row, 'POS'])+'_'+str(df.loc[row, 'REF'])+'_'+str(df.loc[row, 'ALT']))
        row += 1 
    return MID_lst

def COMMON_MEMBERS(a, b):    
    a_set = set(a)
    b_set = set(b)         
    return(list(a_set.intersection(b_set)))  
    

def OBJCHASSIS(keys, chassis):
    obj = {}
    for i in keys:
        other = copy.copy(chassis)
        obj[i] = other
    return obj

def ATTRIBUTE_BASES(chassis, df_path, chromosome, reference, MIDs_lst):
    df = pd.read_csv(df_path, sep='\t')
    ind = 0
    while ind < len(df):
        chassis[MIDs_lst[ind]]['ALT'] = df.iloc[ind,3]
        chassis[MIDs_lst[ind]]['REF'] = df.iloc[ind,2]
        chassis[MIDs_lst[ind]]['POS'] = df.iloc[ind,1]
        chassis[MIDs_lst[ind]]['5P'] = reference[chromosome][int(df.iloc[ind,1])-2]
        chassis[MIDs_lst[ind]]['3P'] = reference[chromosome][int(df.iloc[ind,1])+(len(df.iloc[ind,2])-1)]
        chassis[MIDs_lst[ind]]['CONTEXTE'] = reference[chromosome][int(df.iloc[ind,1])-2] + '[' + df.iloc[ind,2] + '>' + df.iloc[ind,3] + ']' + reference[chromosome][int(df.iloc[ind,1])]
        ind += 1
    return chassis

def ATTRIBUTE_BASES_PENTA(chassis, df_path, chromosome, reference, MIDs_lst):
    df = pd.read_csv(df_path, sep='\t')
    ind = 0
    while ind < len(df):
        chassis[MIDs_lst[ind]]['ALT'] = df.iloc[ind,3]
        chassis[MIDs_lst[ind]]['REF'] = df.iloc[ind,2]
        chassis[MIDs_lst[ind]]['POS'] = df.iloc[ind,1]
        chassis[MIDs_lst[ind]]['5P'] = reference[chromosome][int(df.iloc[ind,1])-2]
        chassis[MIDs_lst[ind]]['3P'] = reference[chromosome][int(df.iloc[ind,1])+(len(df.iloc[ind,2])-1)]
        chassis[MIDs_lst[ind]]['CONTEXTE'] = reference[chromosome][int(df.iloc[ind,1])-3] + reference[chromosome][int(df.iloc[ind,1])-2] + '[' + df.iloc[ind,2] + '>' + df.iloc[ind,3] + ']' + reference[chromosome][int(df.iloc[ind,1])] + reference[chromosome][int(df.iloc[ind,1])+1]
        ind += 1
    return chassis


def ATTRIBUTE_NOT_CLUSTERED(frame, df_path, MIDs_lst):
    df = pd.read_csv(df_path, sep='\t')
    ind = 0
    if len(df) > 1:
        while ind < len(df):
            if ind < len(df)-1:
                if abs(df.iloc[ind,1]-df.iloc[abs(ind-1),1]) > 100 and abs(df.iloc[ind,1]-df.iloc[abs(ind+1),1]) > 100:
                    frame[MIDs_lst[ind]]['CLUSTER'] = 'Not_Clustered'
            else:
                if abs(df.iloc[ind,1]-df.iloc[abs(ind-1),1]) > 100:
                    frame[MIDs_lst[ind]]['CLUSTER'] = 'Not_Clustered'
            ind +=1
    else:
        frame[MIDs_lst[ind]]['CLUSTER'] = 'Not_Clustered'
    return frame

def ATTRIBUTE_CLUSTERED(built, df_path, MIDs_lst, chromosome):
    nb = 0
    df = pd.read_csv(df_path, sep='\t')
    ind = 0
    while ind < len(df):
        if ind < len(df)-1:
            if built[MIDs_lst[ind]]['CLUSTER'] == 'null':
                #First mutation of a new cluster
                if abs(df.iloc[ind,1]-df.iloc[abs(ind-1),1]) > 100 and abs(df.iloc[ind,1]-df.iloc[abs(ind+1),1]) <= 100:
                    nb +=1 
                    common = []
                    cluster = chromosome + '-' + str(nb)
                    built[MIDs_lst[ind]]['CLUSTER'] = cluster
                    common.append(MIDs_lst[ind])
                #Mutation in the middle of a cluster
                if abs(df.iloc[ind,1]-df.iloc[abs(ind-1),1]) <= 100 and abs(df.iloc[ind,1]-df.iloc[abs(ind+1),1]) <= 100:
                    #if first mutation of chromosome
                    if ind == 0:
                        cluster = chromosome + '-' + str(1)
                        common = []
                    built[MIDs_lst[ind]]['CLUSTER'] = cluster
                    common.append(MIDs_lst[ind])
                #Last mutation of a cluster
                if abs(df.iloc[ind,1]-df.iloc[abs(ind-1),1]) <= 100 and abs(df.iloc[ind,1]-df.iloc[abs(ind+1),1]) > 100:
                    built[MIDs_lst[ind]]['CLUSTER'] = cluster
                    common.append(MIDs_lst[ind])
                    for mid in common:
                        built[mid]['COMMON'] = common
        #Last mutation of chromosome
        else:
            if built[MIDs_lst[ind]]['CLUSTER'] == 'null':
                built[MIDs_lst[ind]]['CLUSTER'] = cluster
                common.append(MIDs_lst[ind])
                for mid in common:
                    built[mid]['COMMON'] = common
        ind+=1
    return built
    

def CHANGE_INDEX(df, lst):
    '''
    This function changes the index/row names of a dataframe
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<df>>: Dataframe
    <<lst>>: List of the new indexes
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	df -> Dataframe with new indexes
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 7, 2022
	'''
    df = df.set_axis(lst, axis=0)
    return df

def SELECT_ROWS_BASED_ON_INDEX(df, lst):
    '''
    This function selects rows of a dataframe by index/row names
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<df>>: Dataframe
    <<lst>>: List of indexes to select
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	df -> Subset dataframe
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 7, 2022
	'''
    df = df.loc[lst]
    return df

def CLUSTERED_MIDS(obj):
    '''
    This function retrieves the clustered MIDs from cluster objects
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
	<<obj/dict>>: cluster object
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	list -> List of MIDs in clusters
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 7, 2022
	'''
    lst = []
    for i in obj:
        if obj[i]['CLUSTER'] != 'Not_Clustered' and obj[i]['CLUSTER'] != 'Reads_dont_confirm':
            lst.append(i)
    return lst

#The following functions were originally written for the purposes of the:
#########################################################################
#                      object2vcf(solitary).py                          #
#########################################################################

def SOLITARY_MIDS(obj):
    '''
    This function retrieves the solitary MIDs from cluster objects
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
	<<obj/dict>>: cluster object
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	list -> List of solitary MIDs
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 10, 2022
	'''
    lst = []
    for i in obj:
        if obj[i]['CLUSTER'] == 'Not_Clustered' or obj[i]['CLUSTER'] == 'Reads_dont_confirm':
            lst.append(i)
    return lst

#The following functions were originally written for the purposes of the:
#########################################################################
#                        96MATRIX_GROUPING.py                           #
#########################################################################

def COLUMNS2SINGLE_LIST(df, rowcolumn="MutationType"):
    '''
    This function sums all the numeric columns of a pandas DataFrame
    and converts the result to a list. 
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<df>>: Dataframe
    <<rowcolumn>>: Title of the non numeric column 
                   (default: MutationType)
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	lst -> Sums of all lines of DataFrame (List) 
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    october 27, 2022
	'''
    df.drop(columns = rowcolumn, axis = 1, inplace = True) 
    summed = df.sum(axis=1)
    lst = summed.to_list()
    return lst

#The following functions were originally written for the purposes of the:
#########################################################################
#                           96matrixMerge.py                            #
#########################################################################

def POP_FIRST(df):
    '''
    This function removes the first column of a dataframe  
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<df>>: Dataframe 
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	df -> Dataframe with first column removed  
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    november 2, 2022
	'''
    df = df.iloc[:,1:]
    return df



def DATAFRAME_MERGE(pDF, annDF, pop_first = False):
    '''
    This function annexes a DataFrame to another 
    (equivalent to cbind in R) 
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<pDF>>: principal Dataframe
    <<annDF>>: Dataframe to annexe to pDF
    <<pop_first>>: Boolean value to indicate if the first column
                   needs to be removed (default is False) 
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	pDF -> Annexed DataFrames  
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    november 2, 2022
	'''
    if pop_first == True:
        annDF = POP_FIRST(annDF)
    for title in annDF:
        lst = annDF[title].to_list()
        pDF[title] = lst
    return pDF

#The following functions were originally written for the purposes of the:
#########################################################################
#                        COMPENDIUM2PROFILE.py                          #
#########################################################################

def CLUSTER_WITHOUT_BAM(df, dist = 100, chr_pos = 0, pos_pos = 1, ref_pos = 2, alt_pos = 3, prime5 = 4, prime3 = 5):
    contexte = []
    treshold = len(df)
    cpt = 0
    check = []
    while cpt+1 < treshold:
        if df.iloc[cpt, chr_pos] == df.iloc[cpt+1, chr_pos] and df.iloc[cpt+1, pos_pos] - df.iloc[cpt, pos_pos] <= dist:
            id1 = str(df.iloc[cpt, pos_pos]) + '_' + str(df.iloc[cpt, chr_pos])
            id2 = str(df.iloc[cpt+1, pos_pos]) + '_' + str(df.iloc[cpt+1, chr_pos])
            if id1 not in check:
                check.append(id1)
                contexte.append(df.iloc[cpt, prime5] + '['+ df.iloc[cpt, ref_pos] + '>' + df.iloc[cpt, alt_pos] + ']' + df.iloc[cpt, prime3])
            check.append(id2)
            contexte.append(df.iloc[cpt+1, prime5] + '['+ df.iloc[cpt+1, ref_pos] + '>' + df.iloc[cpt+1, alt_pos] + ']' + df.iloc[cpt+1, prime3])
        cpt+=1
    return contexte

def PROFILE_WITHOUT_BAM(df, ref_pos = 2, alt_pos = 3, prime5 = 4, prime3 = 5):
    contexte = []
    treshold = len(df)
    cpt = 0
    while cpt < treshold:
        contexte.append(df.iloc[cpt, prime5] + '['+ df.iloc[cpt, ref_pos] + '>' + df.iloc[cpt, alt_pos] + ']' + df.iloc[cpt, prime3])
        cpt+=1
    return contexte

def SOLITARY_WITHOUT_BAM(df, dist = 100, chr_pos = 0, pos_pos = 1, ref_pos = 2, alt_pos = 3, prime5 = 4, prime3 = 5):
    contexte = []
    treshold = len(df)
    cpt = 0
    check = []
    while cpt+1 < treshold:
        id1 = str(df.iloc[cpt, pos_pos]) + '_' + str(df.iloc[cpt, chr_pos])
        id2 = str(df.iloc[cpt+1, pos_pos]) + '_' + str(df.iloc[cpt+1, chr_pos])
        if df.iloc[cpt, chr_pos] == df.iloc[cpt+1, chr_pos] and df.iloc[cpt+1, pos_pos] - df.iloc[cpt, pos_pos] <= dist:
                check.append(id1)
                check.append(id2)
        else:
            if id1 not in check:
                contexte.append(df.iloc[cpt, prime5] + '['+ df.iloc[cpt, ref_pos] + '>' + df.iloc[cpt, alt_pos] + ']' + df.iloc[cpt, prime3])
            if id2 not in check and cpt+1 == treshold-1:
                contexte.append(df.iloc[cpt+1, prime5] + '['+ df.iloc[cpt+1, ref_pos] + '>' + df.iloc[cpt+1, alt_pos] + ']' + df.iloc[cpt+1, prime3])

        cpt+=1
    return contexte

def POP_LAST(df, n):
    '''
    This function removes the n last column(s) of a dataframe  
    This function requires to import pandas module.
    ----------------------------------------------------------------
    INPUT
	<<df>>: Dataframe
    <<n>>: Number of columns to remove. 
	________________________________________________________________

	----------------------------------------------------------------
	RETURNS
	df -> Dataframe with n last column(s) removed  
	________________________________________________________________
    
	written by Félix Racine-Brassard
	UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    november 8, 2022
	'''
    df = df.iloc[:,:-n]
    return df


def SELECT_ROWS_BY_VALUE_GREATER_THAN(df, column_to_check, target, inclusive=True):
    '''
    This function subsets rows based on the content of specific row
    in a specific column.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<df>>: Dataframe
         <<column_to_check>>: Name of the column to iterate over.
         <<target>>: Value/String used to select rows. 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> interest: DataFrame containing only the wanted 
                               rows.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    if inclusive == True:
        interest = df[df[column_to_check] >= target]
    else:
        interest = df[df[column_to_check] > target]
    return interest



def FIND_ELEMENTS_IN_LIST_BY_INDEX(lst1, lst2):
    '''
    This function searches elements .
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<lst1>>: list to search
         <<lst2>>: list of indexes used to search 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
        List of selected elements
    _________________________________________________________________
    
    written by GeeksforGeeks
    '''
    return [lst1[i] for i in lst2]


def SELECT_ROWS_BY_VALUE_LOWER_THAN(df, column_to_check, target, inclusive=True):
    '''
    This function subsets rows based on the content of specific row
    in a specific column.
    This function does not require to load any modules.
    ----------------------------------------------------------------
    INPUT
         <<df>>: Dataframe
         <<column_to_check>>: Name of the column to iterate over.
         <<target>>: Value/String used to select rows. 
    ________________________________________________________________

    -----------------------------------------------------------------
    RETURNS
         variable -> interest: DataFrame containing only the wanted 
                               rows.
    _________________________________________________________________
    
    written by Félix Racine-Brassard
    UNIV. DE SHERBROOKE, Pr. Pierre-Étienne Jacques Laboratory
    august 18, 2022
    '''
    if inclusive == True:
        interest = df[df[column_to_check] <= target]
    else:
        interest = df[df[column_to_check] < target]
    return interest


def COLUMNS2AVERAGE_LIST(df, cols, row_index):
    average = []
    dff = COLUMNSELECTION(df, cols)
    for i in row_index:
        uniq_row = dff.loc[i,:].values.tolist()
        average.append(np.mean(uniq_row))
    return average
        

def IS_SBS(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 1
    else:
        return 0

def IS_DBS(ref, alt):
    if len(ref) == 2 and len(alt) == 2:
        return 10
    else:
        return 0

def IS_INDEL(ref, alt):
    if len(ref) != len(alt):
        return 100
    else:
        return 0

def IS_LBS_LONGMUTATION(ref, alt):
    if len(ref) == len(alt) and len(ref)>2:
        return 1000
    else:
        return 0

def SBS_COUNT(obj):
    x = 0
    for i in obj:
        if len(obj[i]['REF']) == 1 and len(obj[i]['ALT']) == 1:
            x += 1
    return x

def DBS_COUNT(obj):
    x = 0
    for i in obj:
        if len(obj[i]['REF']) == 2 and len(obj[i]['ALT']) == 2:
            x += 1
    return x

def INDEL_COUNT(obj):
    x = 0
    for i in obj:
        if len(obj[i]['REF']) != len(obj[i]['ALT']):
            x += 1
    return x

def INSERTION_COUNT(obj):
    x = 0
    for i in obj:
        if len(obj[i]['REF']) < len(obj[i]['ALT']):
            x += 1
    return x

def DELETION_COUNT(obj):
    x = 0
    for i in obj:
        if len(obj[i]['REF']) > len(obj[i]['ALT']):
            x += 1
    return x

def CLUSTERED_COUNT(obj):
    x = 0
    for i in obj:
        if obj[i]['CLUSTER'] != 'Not_Clustered' and obj[i]['CLUSTER'] != 'Reads_dont_confirm':
            x += 1
    return x

def CLUSTER_ANNOTATION(obj):
    tmp = 'null'
    info = {}
    asso_type=[]
    x = set()
    intermut_dist = []
    for i in obj:
        if obj[i]['CLUSTER'] != 'Not_Clustered' and obj[i]['CLUSTER'] != 'Reads_dont_confirm':
            cluster_found = obj[i]['CLUSTER']
            if cluster_found != tmp:
                if len(info) != 0:
                    info[tmp] = {}
                    info[tmp]['NUM'] = n
                    info[tmp]['ASSO'] = asso_type
                    info[tmp]['INTERMUT_DIST'] = intermut_dist
                else:
                    info[tmp] = {} 
                asso_type = []
                intermut_dist = []
                tmp = copy.copy(cluster_found)
                n = 1
            else:
                n += 1
                asso_type.append(IS_SBS(obj[i]['REF'], obj[i]['ALT']) + IS_SBS(obj[itmp]['REF'], obj[itmp]['ALT']) + IS_DBS(obj[i]['REF'], obj[i]['ALT']) + IS_DBS(obj[itmp]['REF'], obj[itmp]['ALT']) + IS_INDEL(obj[i]['REF'], obj[i]['ALT'])+ IS_INDEL(obj[itmp]['REF'], obj[itmp]['ALT']) + IS_LBS_LONGMUTATION(obj[i]['REF'], obj[i]['ALT']) + IS_LBS_LONGMUTATION(obj[itmp]['REF'], obj[itmp]['ALT']))
                intermut_dist.append(obj[itmp]['POS']-(obj[i]['POS']+(len(obj[i]['POS'])-1)))   
            x.add(obj[i]['CLUSTER'])
            itmp = copy.copy(i)
    if len(info) != 0:
        info[tmp] = {}
        info[tmp]['NUM'] = n
        info[tmp]['ASSO'] = asso_type
        info[tmp]['INTERMUT_DIST'] = intermut_dist
        del info['null']
    return x, info

def CLUSTER_PAIRS_COUNT(asso):
    n = 0
    for i in asso:
        if asso[i]['NUM'] == 2:
            n += 1 
    return n

def CLUSTER_TRIPLETS_COUNT(asso):
    n = 0
    for i in asso:
        if asso[i]['NUM'] == 3:
            n += 1 
    return n

def CLUSTER_QUADRUPLETS_COUNT(asso):
    n = 0
    for i in asso:
        if asso[i]['NUM'] == 4:
            n += 1 
    return n

def CLUSTER_FIVE_PLUS_COUNT(asso):
    n = 0
    for i in asso:
        if asso[i]['NUM'] >= 5:
            n += 1 
    return n

def INDEL_SBS_ASSOCIATION_COUNT(asso):
    n = 0
    for i in asso:
        for j in asso[i]['ASSO']:
            if j == 101:
                n += 1 
    return n

def INDEL_INDEL_ASSOCIATION_COUNT(asso):
    n = 0
    for i in asso:
        for j in asso[i]['ASSO']:
            if j == 200:
                n += 1 
    return n

def DBS_DBS_ASSOCIATION_COUNT(asso):
    n = 0
    for i in asso:
        for j in asso[i]['ASSO']:
            if j == 20:
                n += 1 
    return n

def SBS_SBS_ASSOCIATION_COUNT(asso):
    n = 0
    for i in asso:
        for j in asso[i]['ASSO']:
            if j == 2:
                n += 1 
    return n

def DBS_SBS_ASSOCIATION_COUNT(asso):
    n = 0
    for i in asso:
        for j in asso[i]['ASSO']:
            if j == 11:
                n += 1 
    return n

def DISTANCES_SBS_SBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 2:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_SBS_DBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 11:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_SBS_INDEL(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 101:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_SBS_LBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 1001:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_DBS_DBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 20:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_DBS_INDEL(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 110:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_DBS_LBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 1010:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_INDEL_INDEL(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 200:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_INDEL_LBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 1100:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_LBS_LBS(annotation):
    distances = []
    for cluster in annotation:
        for idx, association in enumerate(annotation[cluster]['ASSO']):
            if association == 2000:
                distances.append(annotation[cluster]['INTERMUT_DIST'][idx])
    return distances

def DISTANCES_WITH_DBS_BUT_WITHOUT_INDEL_OR_LBS(annotation):
    distances = DISTANCES_SBS_DBS(annotation) + DISTANCES_DBS_DBS(annotation)
    return distances

def VERIFY_SEQ(sequence):
    '''This code verifies if a sequence is a DNA or RNA'''
    sequp = sequence.upper()
    valid = True
    for char in sequp:
        if char not in ['A', 'T', 'G', 'C', 'U']:
            valid = False
            break
    if valid == False:
        return "Invalid sequence"
    elif 'U' in sequp:
        return "RNA"
    else:
        return "DNA"
    

def REVERSE_COMPLEMENT(seq):
    '''This function returns a reverse complement
    of a DNA or RNA strand'''
    verified = VERIFY_SEQ(seq)
    if verified == "DNA":
       
        # complement strand
        seq = seq.replace("A", "t").replace(
            "C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
         
        # reverse strand
        seq = seq[::-1]
        return seq
 
    elif verified == "RNA":
       
        # complement strand
        seq = seq.replace("A", "u").replace(
            "C", "g").replace("U", "a").replace("G", "c")
        seq = seq.upper()
         
        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return "Invalid sequence"

def CHECK_CONTEXTE_VALIDITY(ctx):
    if '-' in str(ctx) or ctx == 0 or len(str(ctx))>7:
        return False
    else:
        return True

def CHECK_IF_DBS(ref, alt):
    if len(ref) == 2 and len(alt)==2:
        return True
    else:
        return False


def RESOLVE_COMPLEMENT(ctx):
    ctxt96 = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']
    ctxt_comp = {'T[G>T]T': 'A[C>A]A', 'G[G>T]T': 'A[C>A]C', 'C[G>T]T': 'A[C>A]G', 'A[G>T]T': 'A[C>A]T', 'T[G>C]T': 'A[C>G]A', 'G[G>C]T': 'A[C>G]C', 'C[G>C]T': 'A[C>G]G', 'A[G>C]T': 'A[C>G]T', 'T[G>A]T': 'A[C>T]A', 'G[G>A]T': 'A[C>T]C', 'C[G>A]T': 'A[C>T]G', 'A[G>A]T': 'A[C>T]T', 'T[A>T]T': 'A[T>A]A', 'G[A>T]T': 'A[T>A]C', 'C[A>T]T': 'A[T>A]G', 'A[A>T]T': 'A[T>A]T', 'T[A>G]T': 'A[T>C]A', 'G[A>G]T': 'A[T>C]C', 'C[A>G]T': 'A[T>C]G', 'A[A>G]T': 'A[T>C]T', 'T[A>C]T': 'A[T>G]A', 'G[A>C]T': 'A[T>G]C', 'C[A>C]T': 'A[T>G]G', 'A[A>C]T': 'A[T>G]T', 'T[G>T]G': 'C[C>A]A', 'G[G>T]G': 'C[C>A]C', 'C[G>T]G': 'C[C>A]G', 'A[G>T]G': 'C[C>A]T', 'T[G>C]G': 'C[C>G]A', 'G[G>C]G': 'C[C>G]C', 'C[G>C]G': 'C[C>G]G', 'A[G>C]G': 'C[C>G]T', 'T[G>A]G': 'C[C>T]A', 'G[G>A]G': 'C[C>T]C', 'C[G>A]G': 'C[C>T]G', 'A[G>A]G': 'C[C>T]T', 'T[A>T]G': 'C[T>A]A', 'G[A>T]G': 'C[T>A]C', 'C[A>T]G': 'C[T>A]G', 'A[A>T]G': 'C[T>A]T', 'T[A>G]G': 'C[T>C]A', 'G[A>G]G': 'C[T>C]C', 'C[A>G]G': 'C[T>C]G', 'A[A>G]G': 'C[T>C]T', 'T[A>C]G': 'C[T>G]A', 'G[A>C]G': 'C[T>G]C', 'C[A>C]G': 'C[T>G]G', 'A[A>C]G': 'C[T>G]T', 'T[G>T]C': 'G[C>A]A', 'G[G>T]C': 'G[C>A]C', 'C[G>T]C': 'G[C>A]G', 'A[G>T]C': 'G[C>A]T', 'T[G>C]C': 'G[C>G]A', 'G[G>C]C': 'G[C>G]C', 'C[G>C]C': 'G[C>G]G', 'A[G>C]C': 'G[C>G]T', 'T[G>A]C': 'G[C>T]A', 'G[G>A]C': 'G[C>T]C', 'C[G>A]C': 'G[C>T]G', 'A[G>A]C': 'G[C>T]T', 'T[A>T]C': 'G[T>A]A', 'G[A>T]C': 'G[T>A]C', 'C[A>T]C': 'G[T>A]G', 'A[A>T]C': 'G[T>A]T', 'T[A>G]C': 'G[T>C]A', 'G[A>G]C': 'G[T>C]C', 'C[A>G]C': 'G[T>C]G', 'A[A>G]C': 'G[T>C]T', 'T[A>C]C': 'G[T>G]A', 'G[A>C]C': 'G[T>G]C', 'C[A>C]C': 'G[T>G]G', 'A[A>C]C': 'G[T>G]T', 'T[G>T]A': 'T[C>A]A', 'G[G>T]A': 'T[C>A]C', 'C[G>T]A': 'T[C>A]G', 'A[G>T]A': 'T[C>A]T', 'T[G>C]A': 'T[C>G]A', 'G[G>C]A': 'T[C>G]C', 'C[G>C]A': 'T[C>G]G', 'A[G>C]A': 'T[C>G]T', 'T[G>A]A': 'T[C>T]A', 'G[G>A]A': 'T[C>T]C', 'C[G>A]A': 'T[C>T]G', 'A[G>A]A': 'T[C>T]T', 'T[A>T]A': 'T[T>A]A', 'G[A>T]A': 'T[T>A]C', 'C[A>T]A': 'T[T>A]G', 'A[A>T]A': 'T[T>A]T', 'T[A>G]A': 'T[T>C]A', 'G[A>G]A': 'T[T>C]C', 'C[A>G]A': 'T[T>C]G', 'A[A>G]A': 'T[T>C]T', 'T[A>C]A': 'T[T>G]A', 'G[A>C]A': 'T[T>G]C', 'C[A>C]A': 'T[T>G]G', 'A[A>C]A': 'T[T>G]T'}
    if ctx not in ctxt96:
        pyrimidine_ctx = ctxt_comp[ctx]
    else:
        pyrimidine_ctx = ctx
    return pyrimidine_ctx

def RESOLVE_DBS_COMPLEMENT(ref, alt):
    std = ['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG']
    mut = ref + '>' + alt
    if mut not in std:
        mut = REVERSE_COMPLEMENT(ref) + '>' + REVERSE_COMPLEMENT(alt)
    return mut
    
def OBJECT_COUNT_SBS_CONTEXTE(obj, dict_order, type='all'):
    order = copy.copy(dict_order)
    if type == 'all':
        for i in obj:
            if CHECK_CONTEXTE_VALIDITY(obj[i]['CONTEXTE']):
                ctxt = RESOLVE_COMPLEMENT(obj[i]['CONTEXTE'])
                order[ctxt] += 1
    elif type == 'clustered':
        for i in obj:
            if obj[i]['CLUSTER'] != 'Not_Clustered' and obj[i]['CLUSTER'] != 'Reads_dont_confirm':
                if CHECK_CONTEXTE_VALIDITY(obj[i]['CONTEXTE']):
                    ctxt = RESOLVE_COMPLEMENT(obj[i]['CONTEXTE'])
                    order[ctxt] += 1
    elif type == 'solitary':
        for i in obj:
            if obj[i]['CLUSTER'] == 'Not_Clustered' or obj[i]['CLUSTER'] == 'Reads_dont_confirm':
                if CHECK_CONTEXTE_VALIDITY(obj[i]['CONTEXTE']):
                    ctxt = RESOLVE_COMPLEMENT(obj[i]['CONTEXTE'])
                    order[ctxt] += 1
    return order


def MAKE_CONTEXTE_TABLE(contexte_dict_list, mut_order, sample_names):
    df = pd.DataFrame()
    n = 0
    for i in contexte_dict_list:
        lst = []
        for j in mut_order:
            lst.append(i[j])
        df[sample_names[n]] = lst
        n += 1
    df.index = mut_order
    return df

    
def OBJECT_COUNT_DBS_CONTEXTE(obj, dict_order, type='all'):
    order = copy.copy(dict_order)
    if type == 'all':
        for i in obj:
            if CHECK_IF_DBS(obj[i]['REF'], obj[i]['ALT']):
                ctxt = RESOLVE_DBS_COMPLEMENT(obj[i]['REF'], obj[i]['ALT'])
                order[ctxt] += 1
    elif type == 'clustered':
        for i in obj:
            if obj[i]['CLUSTER'] != 'Not_Clustered' and obj[i]['CLUSTER'] != 'Reads_dont_confirm':
                if CHECK_IF_DBS(obj[i]['REF'], obj[i]['ALT']):
                    ctxt = RESOLVE_DBS_COMPLEMENT(obj[i]['REF'], obj[i]['ALT'])
                    order[ctxt] += 1
    elif type == 'solitary':
        for i in obj:
            if obj[i]['CLUSTER'] == 'Not_Clustered' or obj[i]['CLUSTER'] == 'Reads_dont_confirm':
                if CHECK_IF_DBS(obj[i]['REF'], obj[i]['ALT']):
                    ctxt = RESOLVE_DBS_COMPLEMENT(obj[i]['REF'], obj[i]['ALT'])
                    order[ctxt] += 1
    return order

def VCF2tables(vcf_file_path, namespace = "BSgenome.Hsapiens.1000genomes.hs37d5", ref_genome='GRCh37', variant_caller='Mutect2', num_of_cores=1, trans_ranges="trans.ranges.GRCh37", region="genome"):
    base = importr('base')
    icams = importr("ICAMS")
    base.loadNamespace(namespace)
    result = icams.VCFsToCatalogs(StrVector(vcf_file_path), ref_genome, variant_caller,num_of_cores,robjects.r(trans_ranges), region)
    sbsr_df = base.as_data_frame(result.rx2("catSBS96"))
    sbspd_dt = robjects.conversion.rpy2py(sbsr_df)
    sbsdf = pd.DataFrame(sbspd_dt, columns=base.rownames(sbspd_dt), index=base.colnames(sbspd_dt)).T
    dbsr_df = base.as_data_frame(result.rx2("catDBS78"))
    dbspd_dt = robjects.conversion.rpy2py(dbsr_df)
    dbsdf = pd.DataFrame(dbspd_dt, columns=base.rownames(dbspd_dt), index=base.colnames(dbspd_dt)).T
    IDr_df = base.as_data_frame(result.rx2("catID78"))
    IDpd_dt = robjects.conversion.rpy2py(IDr_df)
    IDdf = pd.DataFrame(IDpd_dt, columns=base.rownames(IDpd_dt), index=base.colnames(IDpd_dt)).T
    return sbsdf, dbsdf, IDdf

def SORT_VCF_TABLE(table, chr=False):
    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
    if chr == True:
        chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrMT']
    sorted_table = pd.DataFrame()
    for i in chromosomes:
        tmp = table[table['#CHROM'] == i]
        if tmp.empty == False:
            tmp.sort_values('POS')
        sorted_table = pd.concat([sorted_table, tmp])
    return sorted_table

def REMOVE_VCF_HEADER(vcf_file, store_header = True, chr=False):
    with open(vcf_file, 'r') as file:
        if store_header == True:
            header_lines = []
            colnames = []
            dfrows = []
            for line in file:
                if line.startswith('##'):
                    header_lines.append(line)
                elif line.startswith('#'):
                    title_col = line
                    colnames = line.strip().split('\t')
                else:
                    dfrows.append(line.strip().split('\t'))
            df = pd.DataFrame(dfrows, columns=colnames )
            df = SORT_VCF_TABLE(df, chr=chr)
            return df, title_col, header_lines
        else:
            header_lines = []
            colnames = []
            dfrows = []
            for line in file:
                if line.startswith('##'):
                    pass
                elif line.startswith('#'):
                    colnames = line.strip().split('\t')
                else:
                    dfrows.append(line.strip().split('\t'))
            df = pd.DataFrame(dfrows, columns=colnames )
            return df
        
def RETRIEVE_VAF_FROM_VCF_LINE(row, format='PCAWG'):
    if format == 'PCAWG':
        RowToStr = '\t'.join(map(str, row.to_list()))
        if 'VAF=' in  RowToStr:
            return float(RowToStr[RowToStr.find('VAF=')+4: RowToStr[RowToStr.find('VAF=')+4:].find(';')+RowToStr.find('VAF=')+4]), RowToStr
        else:
            return None, None
                

def OBJ2ISOLATED_CONTEXTS(objet, contexts, objectpassed = False):
    '''.bash_aliases'''
    recipient = {}
    ctx = []
    mid = []
    if objectpassed == False:
        with open(objet, 'rb') as fobj:
            obj = pickle.load(fobj)
    else:
        obj = objet
    for i in obj:
        if any(sub in obj[i]['CONTEXTE'] for sub in contexts): 
            recipient[i] = obj[i]
            ctx.append(obj[i]['CONTEXTE'])
            mid.append(i)
    return recipient, ctx, mid


def GENERATE_CONTEXT_POSSIBILITIES(all_mut_types, which_type='sbs', channel_expansion='false', style='positive_strand', custom='false', context_size = 3):
    '''.bash_aliases'''
    convert_key = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    chnls = []
    r_chnls = []
    chnls_conversion = {}
    flanks = [('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'T'), ('C', 'A'), ('C', 'C'), ('C', 'G'), ('C', 'T'), ('G', 'A'), ('G', 'C'), ('G', 'G'), ('G', 'T'), ('T', 'A'), ('T', 'C'), ('T', 'G'), ('T', 'T')]
    if custom == 'add_flanks':
        for muttype in all_mut_types:
            for cmb in flanks:
                chnls.append(cmb[0] + '[' + muttype + ']' + cmb[1])
                r_muttype = ''
                for base in muttype:
                    if base == '>':
                        r_muttype = r_muttype + base
                    else:
                        r_muttype = r_muttype + convert_key[base]
                r_chnls.append(convert_key[cmb[1]] + '[' + r_muttype + ']' + convert_key[cmb[0]])        
        chnls_conversion = dict(zip(r_chnls, chnls))
        return chnls, chnls_conversion
    if custom == 'complete_contexte':
        if '[' not in all_mut_types and ']':
            return None
    if which_type == 'sbs' and style == 'positive_strand':
        mains = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    elif which_type == 'sbs' and style == 'negative_strand':
        mains = ['G>T', 'G>C', 'G>A', 'A>T', 'A>G', 'A>C'] 
    elif which_type == 'sbs' and channel_expansion == 'True':
        mains = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'G>T', 'G>C', 'G>A', 'A>T', 'A>G', 'A>C']
    elif which_type == 'dbs' and style == 'positive_strand':
        mains = ['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG']
    elif which_type == 'dbs' and style == 'negative_strand':
        print('WARNING: channels AT>CG, AT>GC, AT>TA, CG>AT, CG>GC, CG>TA, GC>AT, GC>CG, GC>TA, TA>AT, TA>CG and TA>GC cannot be differentiated from the positive strand. Hence these channels show the combined count of the positive and negative strand')
        mains = ['GT>TG', 'GT>CG', 'GT>AG', 'GT>TC', 'GT>CC', 'GT>AC', 'GT>TA', 'GT>CA', 'GT>AA', 'AT>TG', 'AT>GG', 'AT>CG', 'AT>TC', 'AT>GC', 'AT>TA', 'GG>TT', 'GG>CT', 'GG>AT', 'GG>TC', 'GG>CC', 'GG>AC', 'GG>TA', 'GG>CA', 'GG>AA', 'CG>AT', 'CG>GC', 'CG>AC', 'CG>TA', 'CG>GA', 'CG>AA', 'AG>TT', 'AG>GT', 'AG>CT', 'AG>TC', 'AG>GC', 'AG>AA', 'AG>TA', 'AG>GA', 'AG>CA', 'GC>TT', 'GC>CT', 'GC>AT', 'GC>TG', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>TG', 'TA>GC', 'TA>CC', 'TA>AC', 'GA>TT', 'GA>CT', 'GA>AT', 'GA>TG', 'GA>CG', 'GA>AG', 'GA>TC', 'GA>CC', 'GA>AC', 'CA>TT', 'CA>GT', 'CA>AT', 'CA>TG', 'CA>GG', 'CA>AG', 'CA>TC', 'CA>GC', 'CA>AC', 'AA>TT', 'AA>GT', 'AA>CT', 'AA>TG', 'AA>GG', 'AA>CG', 'AA>TC', 'AA>GC', 'AA>CC']
    elif which_type == 'dbs' and channel_expansion == 'True':
        print('WARNING: channels AT>CG, AT>GC, AT>TA, CG>AT, CG>GC, CG>TA, GC>AT, GC>CG, GC>TA, TA>AT, TA>CG and TA>GC cannot be differentiated from the positive strand. Hence these channels show the combined count of the positive and negative strand')
        mains = ['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>GA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>GT', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>CA', 'TA>CG', 'TA>CT', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG', 'GT>TG', 'GT>CG', 'GT>AG', 'GT>TC', 'GT>CC', 'GT>AC', 'GT>TA', 'GT>CA', 'GT>AA', 'AT>TG', 'AT>GG', 'AT>CG', 'AT>TC', 'AT>GC', 'AT>TA', 'GG>TT', 'GG>CT', 'GG>AT', 'GG>TC', 'GG>CC', 'GG>AC', 'GG>TA', 'GG>CA', 'GG>AA', 'CG>AT', 'CG>GC', 'CG>AC', 'CG>TA', 'CG>GA', 'CG>AA', 'AG>TT', 'AG>GT', 'AG>CT', 'AG>TC', 'AG>GC', 'AG>AA', 'AG>TA', 'AG>GA', 'AG>CA', 'GC>TT', 'GC>CT', 'GC>AT', 'GC>TG', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>TG', 'TA>GC', 'TA>CC', 'TA>AC', 'GA>TT', 'GA>CT', 'GA>AT', 'GA>TG', 'GA>CG', 'GA>AG', 'GA>TC', 'GA>CC', 'GA>AC', 'CA>TT', 'CA>GT', 'CA>AT', 'CA>TG', 'CA>GG', 'CA>AG', 'CA>TC', 'CA>GC', 'CA>AC', 'AA>TT', 'AA>GT', 'AA>CT', 'AA>TG', 'AA>GG', 'AA>CG', 'AA>TC', 'AA>GC', 'AA>CC']
    elif which_type == 'indel':
        return ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5', '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5', '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5', '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5', '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5', '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5', '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5', '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3', '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5']    
    for muttype in mains:
        for cmb in flanks:
            chnls.append(cmb[0] + '[' + muttype + ']' + cmb[1])
            r_muttype = ''
            for base in reversed(muttype):
                if base == '>':
                    r_muttype = r_muttype + base
                else:
                    r_muttype = r_muttype + convert_key[base]
            r_chnls.append(convert_key[cmb[1]] + '[' + r_muttype + ']' + convert_key[cmb[0]])
    chnls_conversion = dict(zip(r_chnls, chnls))    
    return chnls, chnls_conversion


def BRASSBEDPE_TABLE_FORMAT(bedpe_paths):
    all_bedpe = {}
    for bedpe in bedpe_paths:
        df = pd.read_csv(bedpe, sep='\t', header=9)
        df.rename(columns={'# chr1': 'chrom1', 'chr2':'chrom2'}, inplace=True)
        df = df[df['assembly_score'] != '_']
        all_bedpe[bedpe] = df
    return all_bedpe

def BRASSBEDPE_PLOT_FORMAT(all_bedpe):
    samples = []
    categories = []
    types = []
    counts = []
    template_count = {'Deletion 0-10K':0, 'Deletion 10K-1M':0, 'Deletion > 1M':0, 'Inversion 0-10K':0, 'Inversion 10K-1M':0, 'Inversion > 1M':0, 'Tandem-dup. 0-10K':0, 'Tandem-dup. 10K-1M':0, 'Tandem-dup. > 1M':0, 'Translocation':0}
    for bedpe in all_bedpe:
        cat_count = copy.copy(template_count)
        for index, row in all_bedpe[bedpe].iterrows():
            cat = row['svclass']
            size = int(row['bkdist']) 
            if cat == 'deletion' and size >= 0 and size < 10000:
                cat_count['Deletion 0-10K'] += 1
            elif cat == 'deletion' and size >= 10000 and size < 1000000:
                cat_count['Deletion 10K-1M'] += 1
            elif cat == 'deletion' and size >= 1000000:
                cat_count['Deletion > 1M'] += 1
            elif cat == 'inversion' and size >= 0 and size < 10000:
                cat_count['Inversion 0-10K'] += 1
            elif cat == 'inversion' and size >= 10000 and size < 1000000:
                cat_count['Inversion 10K-1M'] += 1
            elif cat == 'inversion' and size >= 1000000:
                cat_count['Inversion > 1M'] += 1  
            elif cat == 'tandem-duplication' and size >= 0 and size < 10000:
                cat_count['Tandem-dup. 0-10K'] += 1
            elif cat == 'tandem-duplication' and size >= 10000 and size < 1000000:
                cat_count['Tandem-dup. 10K-1M'] += 1
            elif cat == 'tandem-duplication' and size >= 1000000:
                cat_count['Tandem-dup. > 1M'] += 1 
            elif cat == 'translocation':
                cat_count['Translocation'] += 1 
        for cat in cat_count:
            samples.append(bedpe)
            categories.append(cat)
            n = 0
            if 'Translocation' != cat:
                types.append(cat[:cat.find(' ')+n])
            else:
                types.append('Translocation')
            counts.append(cat_count[cat])
    tabling = {'Sample':samples, 'SV_CLASS': types, 'SV_Categories': categories, 'Counts':counts}
    table = pd.DataFrame(tabling)
    return table
            
            
            
                     
        
        