import yaml
import itertools


ext = '.BRASS.readnames.txt'
bamext = '.bam'
obamext = '.BRASS.picardFILTERED.bam'


with open('BRASS.assembledReadNames.catalog.yaml') as yamlfile:
    repo = yaml.safe_load(yamlfile)

with open('BamFiltering.Picard.sh', 'w') as outfile:
    outfile.write('ml picard\n')

for sample in repo:
    s = sample[:sample.find('.filtered')]
    read_list = f'{s}{ext}'
    inputbam = f'{s}{bamext}'
    outputbam = f'{s}{obamext}'
    with open(read_list, 'w') as rlout:
        all_readnames = list(itertools.chain.from_iterable(repo[sample].values()))
        for readname in all_readnames:
            rlout.write(f"{readname}\n")
    with open('BamFiltering.Picard.sh', 'a') as outfile:
        outfile.write(f'java -jar $EBROOTPICARD/picard.jar FilterSamReads I={inputbam} O={outputbam} READ_LIST_FILE={read_list} FILTER=includeReadList\n')

    