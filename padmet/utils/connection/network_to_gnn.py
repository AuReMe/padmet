# -*- coding: utf-8 -*-
    """TODO"""
import docopt
import os
import sys

from cobra.io import read_sbml_model
from collections import OrderedDict
from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def network_to_gnn_cli(command_args):
    """TODO"""


def create_graph_from_padmet(input_file, verbose=False):
    """ Create JSON formatted for metexploreviz using a padmet file

    Parameters
    ----------
    input_file: str
        path to padmet input file
    verbose: bool
        if True print information

    Returns
    -------
    json_dicts: dict
        JSON formatted for metexploreviz
    """
    nodes_in_graph = OrderedDict()
    #loading padmetSpec
    if verbose:
        print('Loading %s' %input_file)
    padmetSpec = PadmetSpec(input_file)

    edges_data = [[], []]
    for node in padmetSpec.dicOfNode.values():
        if node.type == 'reaction':
            reaction_id = node.id
            if node.misc["DIRECTION"][0] == 'REVERSIBLE':
                reversibility = True
            else:
                reversibility = False
            if node.id not in nodes_in_graph:
                nodes_in_graph[reaction_id] = len(nodes_in_graph)
            for rlt in padmetSpec.dicOfRelationIn[reaction_id]:
                if rlt.type in ['consumes']:
                    reactant_id = rlt.id_out
                    if reactant_id not in nodes_in_graph:
                        nodes_in_graph[reactant_id] = len(nodes_in_graph)
                    edges_data[0].append(nodes_in_graph[reactant_id])
                    edges_data[1].append(nodes_in_graph[reaction_id])
                    if reversibility is True:
                        edges_data[1].append(nodes_in_graph[reactant_id])
                        edges_data[0].append(nodes_in_graph[reaction_id])

                if rlt.type in ['produces']:
                    product_id = rlt.id_out
                    if product_id not in nodes_in_graph:
                        nodes_in_graph[product_id] = len(nodes_in_graph)
                    edges_data[0].append(nodes_in_graph[reaction_id])
                    edges_data[1].append(nodes_in_graph[product_id])
                    if reversibility is True:
                        edges_data[1].append(nodes_in_graph[reaction_id])
                        edges_data[0].append(nodes_in_graph[product_id])

    return nodes_in_graph, edges_data


def create_graph_from_sbml(input_file, verbose=False):
    """ Create JSON formatted for metexploreviz using a sbml file

    Parameters
    ----------
    input_file: str
        path to sbml input file
    verbose: bool
        if True print information

    Returns
    -------
    json_dicts: dict
        JSON formatted for metexploreviz
    """
    nodes_in_graph = OrderedDict()

    if verbose:
        print('Loading %s' %input_file)
    sbml_model = read_sbml_model(input_file)

    edges_data = [[], []]
    for rxn in sbml_model.reactions:
        reaction_id = rxn.id
        reversibility = rxn.reversibility

        if reaction_id not in nodes_in_graph:
            nodes_in_graph[reaction_id] = len(nodes_in_graph)
        for reactant in rxn.reactants:
            reactant_id = reactant.id
            if reactant_id not in nodes_in_graph:
                nodes_in_graph[reactant_id] = len(nodes_in_graph)
            edges_data[0].append(nodes_in_graph[reactant_id])
            edges_data[1].append(nodes_in_graph[reaction_id])
            if reversibility is True:
                edges_data[1].append(nodes_in_graph[reactant_id])
                edges_data[0].append(nodes_in_graph[reaction_id])

        for product in rxn.products:
            product_id = product.id
            if product_id not in nodes_in_graph:
                nodes_in_graph[product_id] = len(nodes_in_graph)
            edges_data[0].append(nodes_in_graph[reaction_id])
            edges_data[1].append(nodes_in_graph[product_id])
            if reversibility is True:
                edges_data[1].append(nodes_in_graph[reaction_id])
                edges_data[0].append(nodes_in_graph[product_id])

    return nodes_in_graph, edges_data



def network_to_gnn(input_file_folder, output_file, focus, verbose=False):
    """ Create JSON formatted for metexploreviz using an inptu file or a folder

    Parameters
    ----------
    input_file_folder: str
        path to input file (either padmet or sbml) or a folder containing padmet or sbml
    output_file: str
        path to JSON formatted for metexploreviz
    verbose: bool
        if True print information

    """
    if os.path.isdir(input_file_folder):
        input_type = "dir"
    elif os.path.isfile(input_file_folder):
        input_type = "file"
    else:
        raise TypeError("%s is not a dir or a file or is not accessible." %(input_file_folder))

    _, file_extension = os.path.splitext(input_file_folder)
    json_dicts_tmp = {}
    json_dicts_tmp['nodes'] = []
    json_dicts_tmp['links'] = []
    if file_extension == '.padmet':
        nodes_data, edges_data = create_graph_from_padmet(input_file_folder, verbose)
    elif file_extension == '.sbml':
        nodes_data, edges_data = create_graph_from_sbml(input_file_folder, verbose)
    else:
        print('Incorrect format for {} Need .padmet of .sbml file.'.format(input_file_folder))
        sys.exit()

    return nodes_data, edges_data
