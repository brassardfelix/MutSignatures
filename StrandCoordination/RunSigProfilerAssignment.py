from SigProfilerAssignment import Analyzer as Analyze
import argparse
import sys
import getopt


def main(argv):
    genome = "" 
    inputType = "" 
    outputDir = "" 
    ctx = "" 
    exclude = ""
    probs = "" 
    samples = ""


    try:
        opts, args = getopt.getopt(argv,"hg:i:o:c:e:p:s:")
    except getopt.GetoptError:
        print ('RunSigProfilerAssignment.py -h <HELP> -g <GENOME> -i <INPUT_TYPE> -o <OUTPUT_DIRECTORY> -c <CONTEXT_TYPE> -e <SIG_TYPES_TO_EXCLUDE> -p <BOOLEAN_GET_PROBS> -s <SAMPLES>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('RunSigProfilerAssignment.py -h <HELP> -g <GENOME> -i <INPUT_TYPE> -o <OUTPUT_DIRECTORY> -c <CONTEXT_TYPE> -e <SIG_TYPES_TO_EXCLUDE> -p <BOOLEAN_GET_PROBS> -s <SAMPLES>')
            sys.exit()
        elif opt in ("-g", "--genome"):
            genome = arg
        elif opt in ("-i", "--input-dir"):
            inputType = arg
        elif opt in ("-o", "--output-dir"):
            outputDir = arg
        elif opt in ("-c", "--contextType"):
            ctx = arg
        elif opt in ("-e", "--exclude"):
            exclude = arg
        elif opt in ("-p", "--probabilities"):
            probs = arg
        elif opt in ("-s", "--samples"):
            samples = arg
    
    
    sample_lst = samples.split('%!')
    exclusion = exclude.split('%!')
    Analyze.cosmic_fit(sample_lst, outputDir, input_type=inputType, context_type=ctx,
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                    genome_build=genome, signature_database=None,
                    exclude_signature_subgroups=exclusion, export_probabilities=probs,
                    export_probabilities_per_mutation=probs, make_plots=False,
                    sample_reconstruction_plots=False, verbose=False)
    
if __name__ == '__main__':
    main(sys.argv[1:])
