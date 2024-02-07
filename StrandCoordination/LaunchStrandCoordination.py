from SigProfilerTopography import Topography as topography
import argparse
import sys
import getopt


### Base required parameters
def main(argv):

    genome = "" ## Needs to be ajusted in function of needs
    inputDir = "" ## Needs to be ajusted in function of needs
    outputDir = "" ## Needs to be ajusted in function of needs
    jobname = "" ## Needs to be ajusted in function of needs
    numofSimulations = 5 ## (Default) Can be ajusted in function of needs

    ### Additional parameters ###

    sbs_matrix_path = 'TODO' ## Needs to be ajusted in function of needs
    samples = '' ## Needs to be ajusted in function of needs

    try:
        opts, args = getopt.getopt(argv,"hg:i:o:j:n:m:s",["vcf_file=", "samp="])
    except getopt.GetoptError:
        print ('LaunchStrandCoordination.py -h <HELP> -g <GENOME> -i <INPUT_DIRECTORY> -o <OUTPUT_DIRECTORY> -j <JOB_NAME> -n <NUM_OF_SIMULATIONS> -m <SBS_MATRIX_PATH> -s <SAMPLES>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('LaunchStrandCoordination.py -h <HELP> -g <GENOME> -i <INPUT_DIRECTORY> -o <OUTPUT_DIRECTORY> -j <JOB_NAME> -n <NUM_OF_SIMULATIONS> -m <SBS_MATRIX_PATH> -s <SAMPLES>')
            sys.exit()
        elif opt in ("-g", "--genome"):
            genome = arg
        elif opt in ("-i", "--input-dir"):
            inputDir = arg
        elif opt in ("-o", "--output-dir"):
            outputDir = arg
        elif opt in ("-j", "--job"):
            jobname = arg
        elif opt in ("-n", "--nsims"):
            numofSimulations = arg
        elif opt in ("-m", "--matrix"):
            sbs_matrix_path = arg
        elif opt in ("-s", "--samples"):
            samples = arg
    

    sample_lst = samples.split(';')
    print(f'Strand-Coordination analysis has started for {jobname}.\n\tReference genome is {genome}.\n\tInput directory is {inputDir} and relevant samples are {sample_lst}.\n\tOutput directory is {outputDir}.\n\tAnalyses wil be performed using the signatures shown in {sbs_matrix_path}.\n\t***ALL DBS AND SV COSMIC SIGNATURES ARE USED***\n\t{numofSimulations} simulations will be performed.\n')

    topography.runAnalyses(genome, 
        inputDir, 
        outputDir, 
        jobname, 
        numofSimulations, 
        processivity=True,
        discreet_mode=True, #Default
        plot_processivity = True,
        sbs_signatures = sbs_matrix_path,
        samples_of_interest = sample_lst
        )

if __name__ == '__main__':
    main(sys.argv[1:])