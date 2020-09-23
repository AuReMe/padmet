# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in 1-n or 2 sbml.

    Returns if a reaction is missing

    And if a reaction with the same id is using different species or different reversibility

::

    usage:
        padmet compare_sbml --sbml=FILES/DIR --output=DIR

    option:
        -h --help    Show help.
        --sbml FILES/DIR    pathname of the sbml files, sep all files by ',', ex: /path/sbml1.sbml;/path/sbml2.sbml OR a folder
        --output DIR    pathname of the output folder
"""
import docopt
import csv
import os
import re

from cobra.io.sbml import read_sbml_model
from padmet.utils.sbmlPlugin import convert_from_coded_id
from padmet.utils.gbr import compile_input


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def compare_sbml_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    sbml_path = args["--sbml"]
    output_folder = args["--output"]
    compare_multiple_sbml(sbml_path, output_folder)


def compare_multiple_sbml(sbml_path, output_folder):
    """
    Compare 1-n sbml, create two output files reactions.tsv and metabolites.tsv
    with the reactions/metabolites in each sbml

    Parameters
    ----------
    sbml_path: str
        path to a folder containing sbmls or multiple sbml paths separated by a ','
    output_folder: str
        path to the output folder
    """
    if not os.path.exists(output_folder):
        print("Creating %s" %output_folder)
        os.makedirs(output_folder)
    else:
        print("%s already exist, old comparison output folders will be overwritten" %output_folder)

    if os.path.isdir(sbml_path):
        if not os.path.exists(sbml_path):
            raise FileNotFoundError("No SBML directory (--sbml/sbml_path) accessible at " + sbml_path)
        all_files = [os.path.join(sbml_path, f) for f in next(os.walk(sbml_path))[2]]
    else:
        all_files = sbml_path.split(",")
        for sbml_file in all_files:
            if not os.path.exists(sbml_file):
                raise FileNotFoundError("No SBML file (--sbml/sbml_path) accessible at " + sbml_file)

    species_columns = [os.path.splitext(os.path.basename(all_file))[0] for all_file in sorted(all_files)]
    gene_columns = [os.path.splitext(os.path.basename(all_file))[0] + '_genes_assoc (sep=;)' for all_file in sorted(all_files)]
    all_reactions = {}
    all_compounds = []
    reactions = {}
    compounds = {}
    for sbml_file in all_files:
        sbml_1 = read_sbml_model(sbml_file)
        reactions[sbml_file] = sbml_1.reactions
        for rxn in sbml_1.reactions:
            if rxn.id not in all_reactions:
                all_reactions[rxn.id] = rxn
        compounds[sbml_file] = [metabolite.id for metabolite in sbml_1.metabolites]
        all_compounds.extend([metabolite.id for metabolite in sbml_1.metabolites])

    all_compounds = set(all_compounds)

    reaction_file = output_folder + '/reactions.tsv'
    reaction_file_rows = []

    for reaction_id in all_reactions:
        reaction_presents = []
        reaction_genes = []
        row = [reaction_id]
        for sbml_file in sorted(all_files):
            if reaction_id in [rxn.id for rxn in reactions[sbml_file]]:
                reaction_presents.append(1)
            else:
                reaction_presents.append(0)
            if reaction_id in reactions[sbml_file]:
                species_reaction = reactions[sbml_file].get_by_id(reaction_id)
                if 'GENE_ASSOCIATION' in species_reaction.notes:
                    ga_for_gbr = species_reaction.notes['GENE_ASSOCIATION']
                    ga_for_gbr = re.sub(r" or " , "|", ga_for_gbr)
                    ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
                    ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
                    if re.findall("\||&", ga_for_gbr):
                        to_compare_ga_subsets = list(compile_input(ga_for_gbr))
                        genes = []
                        for to_compare_subset in to_compare_ga_subsets:
                            for gene in to_compare_subset:
                                genes.append(gene)
                    else:
                        genes  = [ga_for_gbr.replace('(', '').replace(')', '')]
                    reaction_genes.append(';'.join(genes))
                else:
                    reaction_genes.append('')
            else:
                reaction_genes.append('')

        row = row + reaction_presents + reaction_genes
        row.append(all_reactions[reaction_id].reaction)

        reaction_file_rows.append(row)

    with open(reaction_file, 'w') as output_reaction:
        csvwriter = csv.writer(output_reaction, delimiter='\t')
        csvwriter.writerow(['reaction', *species_columns, *gene_columns, '_formula'])
        csvwriter.writerows(reaction_file_rows)

    compounds_file = output_folder + '/metabolites.tsv'
    compounds_rows = []

    for compound_id in all_compounds:
        row = [compound_id]
        for sbml_file in sorted(all_files):
            if compound_id in compounds[sbml_file]:
                row.append(1)
            else:
                row.append(0)
        compounds_rows.append(row)

    with open(compounds_file, 'w') as output_compound:
        csvwriter = csv.writer(output_compound, delimiter='\t')
        csvwriter.writerow(['metabolite', *sorted(all_files)])
        csvwriter.writerows(compounds_rows)


def compare_sbml(sbml1_path, sbml2_path):
    """
    Compare 2 sbml, print nb of metabolites and reactions.
    If reaction missing print reaction id, and reaction formula.

    Parameters
    ----------
    sbml1_path: str
        path to the first sbml file to compare
    sbml2_path: str
        path to the second sbml file to compare
    """

    if not os.path.exists(sbml1_path):
        raise FileNotFoundError("No SBML file or directory (sbml1_path) accessible at " + sbml1_path)

    if not os.path.exists(sbml2_path):
        raise FileNotFoundError("No SBML file or directory (sbml2_path) accessible at " + sbml2_path)

    sbml_1 = read_sbml_model(sbml1_path)
    sbml_2 = read_sbml_model(sbml2_path)
    
    print("sbml1:")
    print("metabolites: %s" %(len(sbml_1.metabolites)))
    print("reactions: %s" %(len(sbml_1.reactions)))
    print("sbml2:")
    print("metabolites: %s" %(len(sbml_2.metabolites)))
    print("reactions: %s" %(len(sbml_2.reactions)))
    not_in1 = [i for i in sbml_2.reactions if i not in sbml_1.reactions]
    print("reactions not in sbml1: %s" %len(not_in1))
    for i in not_in1:
        print("\t%s" %i)
    not_in2 = [i for i in sbml_1.reactions if i not in sbml_2.reactions]
    print("reactions not in sbml2: %s" %len(not_in2))
    for j in not_in2:
        print("\t%s" %j)
    
    all_diff = set()
    for rxn1 in sbml_1.reactions:
        rxn_id = rxn1.id
        try:
            rxn2 = sbml_2.reactions.get_by_id(rxn_id)
            same_cpd, same_rev =  compare_rxn(rxn1, rxn2)
            if rxn_id not in all_diff:
                if not same_cpd:
                    print("%s use different species" %rxn_id)
                if not same_rev:
                    print("%s use different reversibility" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass
                
    for rxn2 in sbml_2.reactions:
        rxn_id = rxn2.id
        try:
            rxn1 = sbml_1.reactions.get_by_id(rxn_id)
            same_cpd, same_rev =  compare_rxn(rxn1, rxn2)
            if rxn_id not in all_diff:
                if not same_cpd:
                    print("%s use different species" %rxn_id)
                if not same_rev:
                    print("%s use different reversibility" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass

def compare_rxn(rxn1,rxn2):
    """
    compare two cobra reaction object and return (same_cpd, same_rev)
    same_cpd: bool, if true means same compounds consumed and produced
    same_reve: bool, if true means same direction of reaction (reversible or not)

    Parameters
    ----------
    rxn1: cobra.model.reaction
        reaction as cobra object
    rxn2: cobra.model.reaction
        reaction as cobra object

    Returns
    -------
    tuple:
        (same_cpd (bool), same_rev (bool))
    """
    same_cpd, same_rev = False, False
    rxn1_dict,rxn2_dict = {}, {}
    for k,v in list(rxn1.metabolites.items()):
        rxn1_dict[k.id] = v
    for k,v in list(rxn2.metabolites.items()):
        rxn2_dict[k.id] = v
    if rxn1_dict == rxn2_dict: 
        same_cpd = True
    if rxn1.reversibility == rxn2.reversibility:

       same_rev = True

    return (same_cpd, same_rev)

