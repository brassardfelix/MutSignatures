import pandas as pd
from tqdm import tqdm
import json
import gzip
import io
import os

##############################
#### FUNCTION DEFINITIONS ####
##############################

def MERGE2STRING_MISSENSE(c1, c2, c3, c4):
    return str(c1)+'_'+str(c2)+'_'+str(c3)+'_'+str(c4)

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

def REMOVE_GZVCF_HEADER(vcf_file, store_header = True, chr=False):
    with io.TextIOWrapper(gzip.open(vcf_file, 'r')) as file:
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

#######################
#### NEEDED INPUTS ####
#######################

# Files to test
#patient_vcf = ['PACA-CA.DO231284.SA600968.wgs.20201119.gatk-mutect2.somatic.snv.open-filter.vcf.gz'] # Test Data

# Retrieve all vcf.gz files from current working directory
patient_vcf = []
for file in os.listdir():
    if file.endswith(".vcf.gz"):
            patient_vcf.append(file)

# Load pathogenic mutations
with open('/home/local/USHERBROOKE/racf2402/Bureau/GitHub_link/MutSignatures/Cohort_caracterization/ressource_data/FANC.pathogenic.SNVs.json') as jfile:
    missense_repo = json.load(jfile)

# Storing object
Alterations = {}

# Counter
numofdo = 0

for vcf in tqdm(patient_vcf):
    donor_id = vcf[vcf.find('.DO')+1:vcf.find('.SA')]
    df, title_col, header_lines = REMOVE_GZVCF_HEADER(vcf, chr=True)
    donor_muts = df.apply(lambda row: MERGE2STRING_MISSENSE(row['#CHROM'][3:],row['POS'], row['REF'], row['ALT']),axis=1)
    for gene in missense_repo:
        found = set(missense_repo[gene]) & set(list(donor_muts))
        if len(found) != 0:
            numofdo += 1
            Alterations[donor_id] = []
            for alt in found:
                Alterations[donor_id].append((gene, alt))
    if donor_id not in list(Alterations.keys()):
        Alterations[donor_id] = 'No altertions reported'

with open('REPORTED.FANC.pathogenic.SNVs.json', 'w') as outfile:
    json.dump(Alterations, outfile)

print(f'Found {numofdo} donors with pathogenic FANC alterations')

