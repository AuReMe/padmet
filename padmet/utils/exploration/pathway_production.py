# -*- coding: utf-8 -*-
"""
Description:
    Compare 1-n padmet to show the pathway input/output for them.
    This is similar to the Metabolite table at the following address:
    https://biocyc.org/comp-genomics?type=(PWY-CLASS-DIST+%7cAmino-Acid-Biosynthesis%7c)&orgids=(AURANTIMONAS+BSUB+ECOLI+EREC)
::

    usage:
        padmet pathway_production --padmet=FILES/DIR --output=DIR [--cpu INT] [-v]

    option:
        -h --help    Show help.
        --padmet=FILES/DIR    pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
        --output=DIR    pathname of the output folder
        --cpu INT    number of CPU to use in multiprocessing
"""
import docopt
import csv
import os

from multiprocessing import Pool
from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def pathway_production_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    output = args["--output"]
    verbose = args["-v"]
    padmet_path = args["--padmet"]
    number_cpu = args["--cpu"]
    pathway_production(padmet_path, output, verbose, number_cpu)


def extract_pahways_inputs_outputs(file_path):
    """ Extract pathway input and output.

    Parameters
    ----------
    file_path: str
        Pathname of the padmet file.

    Returns
    -------
    pathway_inputs: dict
        Input metaboltie of pathways as key and pathways associated to the metabolites as values.
    pathways_outputs: dict
        Ouput metaboltie of pathways as key and pathways associated to the metabolites as values.
    """
    pathway_inputs = {}
    pathways_outputs = {}
    padmet = PadmetSpec(file_path)
    padmet_pathways = [node for node in list(padmet.dicOfNode.values()) if node.type == 'pathway']
    for pathway in padmet_pathways:
        pathway_name = pathway.misc['COMMON-NAME'][0]
        dict_input_name = pathway.id
        if 'INPUT-COMPOUNDS' in pathway.misc:
            for input_compound in pathway.misc['INPUT-COMPOUNDS'][0].split(','):
                if input_compound not in pathway_inputs:
                    pathway_inputs[input_compound] = [dict_input_name]
                else:
                    pathway_inputs[input_compound].append(dict_input_name)
        if 'OUTPUT-COMPOUNDS' in pathway.misc:
            for output_compound in pathway.misc['OUTPUT-COMPOUNDS'][0].split(','):
                if output_compound not in pathways_outputs:
                    pathways_outputs[output_compound] = [dict_input_name]
                else:
                    pathways_outputs[output_compound].append(dict_input_name)

    return pathway_inputs, pathways_outputs


def pathway_production(padmet_path, output, verbose=None, number_cpu=None):
    """ Create two files degradation_matrix.tsv and biosynthesis_matrix.tsv.
    These files have metabolite as row and organism as column.
    It shows the input (degradation_matrix.tsv) and output (biosynthesis_matrix.tsv) of pathways in the organism.

    Parameters
    ----------
    padmet_path: str
        pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
    output: str
        pathname of the output folder
    verbose: bool
        if True print information
    number_cpu: bool
        Number of CPU
    """
    if not os.path.exists(output):
        if verbose: print("Creating %s" %output)
        os.makedirs(output)
    else:
        if verbose: print("%s already exist, old comparison output folders will be overwritten" %output)

    if os.path.isdir(padmet_path):
        all_files = [os.path.join(padmet_path, f) for f in next(os.walk(padmet_path))[2]]
    else:
        all_files = padmet_path.split(",")

    if number_cpu:
        try:
            number_cpu_to_use = int(number_cpu)
        except ValueError:
            raise ValueError('The number of CPU must be an integer.')
    else:
        number_cpu_to_use = 1

    all_metabolites = []
    all_pathways = {}
    organisms = []
    for padmet_file_path in all_files:
        padmet_id = os.path.splitext(os.path.basename(padmet_file_path))[0]
        pathway_inputs, pathways_outputs = extract_pahways_inputs_outputs(padmet_file_path)
        all_pathways[padmet_id] = pathway_inputs, pathways_outputs
        all_metabolites.extend(pathway_inputs.keys())
        all_metabolites.extend(pathways_outputs.keys())
        organisms.append(padmet_id)

    all_metabolites = set(all_metabolites)

    degradation_matrix = []
    biosynthesis_matrix = []
    for metabolite in all_metabolites:
        degradation_matrix.append([metabolite] + [','.join(all_pathways[organism][0][metabolite]) if metabolite in all_pathways[organism][0] else '' for organism in organisms])
        biosynthesis_matrix.append([metabolite] + [','.join(all_pathways[organism][1][metabolite]) if metabolite in all_pathways[organism][1] else '' for organism in organisms])

    degradation_file = os.path.join(output, 'degradation_matrix.tsv')
    with open(degradation_file, 'w') as degradation_output_file:
        csvwriter = csv.writer(degradation_output_file, delimiter='\t')
        csvwriter.writerow(['Metaboltie', *organisms])
        for row in degradation_matrix:
            if ''.join(row[1:]) != '':
                csvwriter.writerow([*row])

    biosynthesis_file = os.path.join(output, 'biosynthesis_matrix.tsv')
    with open(biosynthesis_file, 'w') as biosynthesis_output_file:
        csvwriter = csv.writer(biosynthesis_output_file, delimiter='\t')
        csvwriter.writerow(['Metaboltie', *organisms])
        for row in biosynthesis_matrix:
            if ''.join(row[1:]) != '':
                csvwriter.writerow([*row])