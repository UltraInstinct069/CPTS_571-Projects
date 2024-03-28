# README

command:
python <python_file> <FASTA_FILE> <0 or 1 for global local alignment respectively> <Config_file>

Example command:

python p1.py Human-Mouse-BRCA2-cds.fasta 0 parameters.config

Output:
Output will be generated in named 'global_output.txt' and 'local_output.txt'

Config file:
parameters.config contains the configuration.

Notes:
For "Human-Mouse-BRCA2-cds", we have selected first 5000 characters.
