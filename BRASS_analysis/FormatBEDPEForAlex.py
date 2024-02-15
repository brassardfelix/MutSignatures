import os
import pandas as pd
os.chdir('/home/racinef/tables/BRASS_BEDPE/')
print(os.getcwd())
final = pd.DataFrame()
files = ['HRN_T6_C12_2-1707693.filtered.bedpe', 'HRN_T6_E4_2-1707694.filtered.bedpe', 'HRN_T6_H12_2-1707695.filtered.bedpe', 'HRU_2-2J_T3_C8_2-1707687.filtered.bedpe', 'HRU_2-2J_T3_G4_2-1707688.filtered.bedpe', 'HRU_UV_2-2J_T3_F12_2-2414933_S1_L002.filtered.bedpe', 'HWN_T6_D1_2-1707692.filtered.bedpe', 'HWN_T6_D7_2-1707691.filtered.bedpe', 'HWN_T6_F2_2-1707690.filtered.bedpe', 'HWU_UV_2-2J_T3_B11_2-2414937_S5_L002.filtered.bedpe', 'HWU_UV_2-2J_T3_D9_2-2414934_S2_L002.filtered.bedpe', 'HWU_UV_2-2J_T3_F11_2-2414935_S3_L002.filtered.bedpe']
for file in files:
    df = pd.read_csv(file, sep='\t')
    df.columns
    tmp = df[['sample','svclass', '# chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'bkdist', 'micro-homology', 'assembled read count', 'occL', 'occH', 'Brass Notation']]
    final = final.append(tmp)

final.to_csv('RFWD3.formatted.bedpe.tsv', sep='\t', index=False)
