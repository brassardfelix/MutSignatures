from CLP_dependency_project_FUNCTIONS import SORT_VCF_TABLE, MERGE_GROUPED_MUTATIONS, REMOVE_VCF_HEADER
import os

vcf_files = os.listdir()
for file in vcf_files:
    df, title, header = REMOVE_VCF_HEADER(file, chr=True)
    dff = MERGE_GROUPED_MUTATIONS(df)
    name = f'{file[:-3]}.mutectFormat.vcf'
    with open(name, 'w') as newfile:
        for line in header:
            newfile.write(line)
        newfile.write(dff.to_csv(index=False, sep='\t'))
            
        
        
    