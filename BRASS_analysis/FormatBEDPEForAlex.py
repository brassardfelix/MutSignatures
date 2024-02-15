import os
import pandas as pd
os.chdir('/home/racinef/tables/BRASS_BEDPE/')
print(os.getcwd())
final = pd.DataFrame()
files = ['HRI_0-07nM_T6_C6_2-1707682.filtered.bedpe', 'HRI_0-07nM_T6_E2_2-1707683.filtered.bedpe', 'HWI_2-1_A6_2-1153454.filtered.bedpe', 'HWI_2-1_B6_2-1153455.filtered.bedpe', 'HWI_2-1_D8_2-1153456.filtered.bedpe']
for file in files:
    df = pd.read_csv(file, sep='\t')
    df.columns
    tmp = df[['sample','svclass', '# chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'bkdist', 'micro-homology', 'assembled read count', 'occL', 'occH', 'Brass Notation']]
    final = final.append(tmp)

final.to_csv('RFWD3.ILS.formatted.bedpe.tsv', sep='\t', index=False)
