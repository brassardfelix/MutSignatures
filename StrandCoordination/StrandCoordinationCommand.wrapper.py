import json

with open('/home/local/USHERBROOKE/racf2402/Bureau/GitHub_link/MutSignatures/StrandCoordination/StrandCoordination.config.json') as config_file:
    parameters = json.load(config_file)

with open('StrandCoordinationLauch.sh', 'w') as file:
    for analysis in parameters:
        samples = ' '.join(parameters[analysis]["samples"])
        file.write(f'python3 LaunchStrandCoordination.py -g {parameters[analysis]["genome"]} -i {parameters[analysis]["InputDir"]} -o {parameters[analysis]["OutputDir"]} -j {parameters[analysis]["jobname"]} -n {parameters[analysis]["nsims"]} -m {parameters[analysis]["matrix"]} -s {samples}\n')