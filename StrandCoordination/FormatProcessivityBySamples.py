import pandas as pd

### test data
### import result file from ip34 (assumes you're in home)
file_location = '/home/racinef/TopoResults/20UV_CtrlSamples/data/processivity/tables/20UV_CtrlSamples_Signatures_Processivity.txt'
signature_processivity = pd.read_csv(file_location, sep='\t', skipfooter=267)
sample_processivity_raw_data = pd.read_csv(file_location, sep='\t', header=10)