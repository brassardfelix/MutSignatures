import pandas as pd
import json

missense_reference_file = '/home/racinef/missense_data/AlphaMissense_hg38.tsv' #Absolute path for ip34
FANC_genes = {'BRCA1':'P38398', 'BRCA2':'P51587', 'ERCC4':'Q92889', 'FANCA':'O15360', 'FANCB':'Q8NB91', 
              'FANCC':'Q00597', 'FANCD2':'Q9BXW9', 'FANCE':'Q9HB96', 'FANCF':'Q9NPI8', 'FANCG':'O15287', 
              'FANCI':'Q9NVI1', 'FANCJ':'Q9BX63', 'FANCL':'Q9NW38', 'FANCM':'Q8IYD8', 'FANCN':'Q86YC2', 
              'FANCW':'Q6PCD5', 'RAD51':'Q06609', 'RAD51C':'O43502', 'REV7':'Q9UI95', 'SLX4':'Q8IY92', 
              'UBE2T':'Q9NPD8', 'XRCC2':'O43543'} #Keys: gene name, Values: Uniprot ID (required for search in alphamissense)

def MERGE2STRING_MISSENSE(c1, c2, c3, c4):
    return str(c1)+'_'+str(c2)+'_'+str(c3)+'_'+str(c4)

missenses = pd.read_csv(missense_reference_file, sep='\t', header=3)

pathogenic_muts = {}
ambiguous_muts = {}
benign_muts = {}

for gene in FANC_genes:
    # Get gene specific missense table
    pathogenics = missenses[(missenses['uniprot_id'] == FANC_genes[gene]) & (missenses['am_class'] == 'likely_pathogenic')]
    ambiguous = missenses[(missenses['uniprot_id'] == FANC_genes[gene]) & (missenses['am_class'] == 'ambiguous')]
    benigns = missenses[(missenses['uniprot_id'] == FANC_genes[gene]) & (missenses['am_class'] == 'ambiguous')]
    # Create mutation IDs
    p_muts = pathogenics.apply(lambda row: MERGE2STRING_MISSENSE(row['#CHROM'][3:],row['POS'], row['REF'], row['ALT']),axis=1)
    a_muts = ambiguous.apply(lambda row: MERGE2STRING_MISSENSE(row['#CHROM'][3:],row['POS'], row['REF'], row['ALT']),axis=1)
    b_muts = benigns.apply(lambda row: MERGE2STRING_MISSENSE(row['#CHROM'][3:],row['POS'], row['REF'], row['ALT']),axis=1)
    # Store mutation IDs
    pathogenic_muts[gene] = list(p_muts)
    ambiguous_muts[gene] = list(a_muts)
    benign_muts[gene] = list(b_muts)

with open('FANC.pathogenic.SNVs.json', 'w') as outfile:
    json.dump(pathogenic_muts, outfile)

with open('FANC.ambiguous.SNVs.json', 'w') as outfile:
    json.dump(ambiguous_muts, outfile)

with open('FANC.benign.SNVs.json', 'w') as outfile:
    json.dump(benign_muts, outfile)