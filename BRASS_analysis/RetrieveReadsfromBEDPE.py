import pandas as pd
import yaml
import os

#############################################################
### PRIOR TO USING THIS SCRIPT, BEDPEs SHOULD BE FILTERED ###
#############################################################

all_dir_files = os.listdir('./')
fbedpe = []
for file in all_dir_files:
    if file.endswith('.filtered.bedpe'):
        fbedpe.append(file)

reads_repo = {}
for bedpe in fbedpe:
    deletion_cnt = 0
    inversion_cnt = 0
    tdup_cnt = 0
    translocation_cnt = 0
    reads_repo[bedpe] = {}
    table = pd.read_csv(bedpe, sep='\t')
    for readnames in table[table['svclass'] == 'deletion']['assembled readnames']:
        deletion_cnt += 1
        reads_repo[bedpe][f'Deletion{deletion_cnt}'] = readnames.split(',')
    for readnames in table[table['svclass'] == 'inversion']['assembled readnames']:
        inversion_cnt += 1
        reads_repo[bedpe][f'Inversion{inversion_cnt}'] = readnames.split(',')
    for readnames in table[table['svclass'] == 'tandem-duplication']['assembled readnames']:
        tdup_cnt += 1
        reads_repo[bedpe][f'Tandem-Duplication{tdup_cnt}'] = readnames.split(',')
    for readnames in table[table['svclass'] == 'tandem-duplication']['assembled readnames']:
        translocation_cnt += 1
        reads_repo[bedpe][f'translocation{translocation_cnt}'] = readnames.split(',')

with open('BRASS.assembledReadNames.catalog.yaml', 'w') as file:
    yaml.dump(reads_repo, file)
            
    



