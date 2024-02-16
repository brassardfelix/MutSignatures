from CLP_dependency_project_FUNCTIONS import REMOVE_VCF_HEADER
from tqdm import tqdm
import sys


def main(files):
    for file in tqdm(files):
        table, cols, header = REMOVE_VCF_HEADER(file)
        name = f'{file[:-4]}.sorted.vcf'
        with open(name, 'w') as newfile:
            for line in header:
                newfile.write(line)
            newfile.write(table.to_csv(index=False, sep='\t'))

if __name__ == '__main__':
    main(sys.argv[1:])