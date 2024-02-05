import pandas as pd
import pysam


### TEST SAMPLE ###
### Assumes that we are located at /home/racinef/projects/def-jacquesp/racinef/BRASS_bedpe of Beluga server
### Import bam file
samfile = pysam.AlignmentFile('../../cancer/HAP1/HRN_T6_C12_2-1707693.bam', "rb" )
### Import bedpe file for sample
bedpe = pd.read_csv('HRN_T6_C12_2-1707693.annot.bedpe', sep='\t') 


