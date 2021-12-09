# -*- coding: utf-8 -*-
"""
Description:
    Compare 1-n padmet to show the pathway input/output for them.
    This is similar to the Metabolite table at the following address:
    https://biocyc.org/comp-genomics?type=(PWY-CLASS-DIST+%7cAmino-Acid-Biosynthesis%7c)&orgids=(AURANTIMONAS+BSUB+ECOLI+EREC)
::

    usage:
        padmet pathway_production --padmet=FILES/DIR --output=DIR [--cpu INT] [--padmetRef=FILE] [-v] [--completion-ratio=FLOAT]

    option:
        -h --help    Show help.
        --padmet=FILES/DIR    pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
        --output=DIR    pathname of the output folder
        --cpu INT    number of CPU to use in multiprocessing
        --padmetRef=FILE    pathanme of the database ref in padmet
        --completion-ratio=FLOAT    float between 0 and 1 used to select pathway with a higher completness
"""
import docopt
import csv
import os
import sys

from multiprocessing import Pool
from padmet.classes import PadmetRef, PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def pathway_production_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    output = args['--output']
    verbose = args['-v']
    padmet_path = args['--padmet']
    padmet_ref_path = args['--padmetRef']
    number_cpu = args['--cpu']
    pathway_completion_ratio = args['--completion-ratio']
    pathway_production(padmet_path, output, verbose, number_cpu, padmet_ref_path, pathway_completion_ratio)


def extract_pahways_inputs_outputs(file_path, padmet_ref_pathways, pathway_completion_ratio):
    """ Extract pathway input and output.

    Parameters
    ----------
    file_path: str
        Pathname of the padmet file.
    padmet_ref_pathways: dict
        Pathways present in padmet ref associated to their reactions.
    pathway_completion_ratio: float
        Ratio given by the user about the completion requires by a pathway to be kept

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

    if padmet_ref_pathways is not None and pathway_completion_ratio is not None:
        nb_rxns_pathways = {pwy.id: set([rlt.id_in for rlt in padmet.dicOfRelationOut.get(pwy.id,[]) if rlt.type == "is_in_pathway"]) for pwy in padmet_pathways}
        ratio_pathways = {pwy_id: len(nb_rxns_pathways[pwy_id])/len(padmet_ref_pathways[pwy_id]) for pwy_id in nb_rxns_pathways if pwy_id in padmet_ref_pathways}
        keep_pathways = {pwy_id for pwy_id in ratio_pathways if ratio_pathways[pwy_id] >= pathway_completion_ratio}
        padmet_pathways = [node for node in padmet_pathways if node.id in keep_pathways]

    for pathway in padmet_pathways:
        pathway_name = pathway.misc['COMMON-NAME'][0]
        pathway_id = pathway.id

        if 'INPUT-COMPOUNDS' in pathway.misc:
            for input_compound in pathway.misc['INPUT-COMPOUNDS'][0].split(','):
                if input_compound not in pathway_inputs:
                    pathway_inputs[input_compound] = [pathway_id]
                else:
                    pathway_inputs[input_compound].append(pathway_id)
        if 'OUTPUT-COMPOUNDS' in pathway.misc:
            for output_compound in pathway.misc['OUTPUT-COMPOUNDS'][0].split(','):
                if output_compound not in pathways_outputs:
                    pathways_outputs[output_compound] = [pathway_id]
                else:
                    pathways_outputs[output_compound].append(pathway_id)

    return pathway_inputs, pathways_outputs


def pathway_production(padmet_path, output, verbose=None, number_cpu=None, padmet_ref_path=None, pathway_completion_ratio=None):
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

    if padmet_ref_path is None and pathway_completion_ratio is not None:
        sys.exit('pathway_completion_ratio option needs a padmetRef to compute the pathway completness ratio.')

    if pathway_completion_ratio is not None:
        try:
            pathway_completion_ratio = float(pathway_completion_ratio)
        except ValueError:
            sys.exit('pathway_completion_ratio must be a float')
        if pathway_completion_ratio < 0 or pathway_completion_ratio > 1:
            sys.exit('pathway_completion_ratio must be < ' + str(0) + ' and > ' + str(1))

    if padmet_ref_path:
        padmet_ref_pathways = {}
        padmetRef = PadmetRef(padmet_ref_path)
        all_pwys = [node for node in list(padmetRef.dicOfNode.values()) if node.type == 'pathway']
        for pwy in all_pwys:
            all_rxns = set([rlt.id_in for rlt in padmetRef.dicOfRelationOut.get(pwy.id,[]) if rlt.type == "is_in_pathway"])
            padmet_ref_pathways[pwy.id] = all_rxns
    else:
        padmet_ref_pathways = None

    all_metabolites = []
    all_pathways = {}
    organisms = []
    for padmet_file_path in all_files:
        padmet_id = os.path.splitext(os.path.basename(padmet_file_path))[0]
        pathway_inputs, pathways_outputs = extract_pahways_inputs_outputs(padmet_file_path, padmet_ref_pathways, pathway_completion_ratio)
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