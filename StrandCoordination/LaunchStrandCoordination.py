from SigProfilerTopography import Topography as topography

### Base required parameters

genome = "GRCh37" ## Needs to be ajusted in function of needs
inputDir = "TODO" ## Needs to be ajusted in function of needs
outputDir = "TODO" ## Needs to be ajusted in function of needs
jobname = "TODO" ## Needs to be ajusted in function of needs
numofSimulations = 5 ## (Default) Can be ajusted in function of needs

### Additional parameters ###

sbs_matrix_path = 'TODO'

topography.runAnalyses(genome, 
    inputDir, 
    outputDir, 
    jobname, 
    numofSimulations, 
    processivity=True,
    discreet_mode=True, #Default
    plot_processivity = True,
    sbs_signatures = sbs_matrix_path
    )