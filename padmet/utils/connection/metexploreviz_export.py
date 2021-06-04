# -*- coding: utf-8 -*-
"""
Description:
    convert a padmet representing a metabolic network into a json compatible with MetExplore.

    usage:
        padmet metexploreviz_export --input=FILE --output=DIR [-v] [--focus=FILE]

    options:
        -h --help     Show help.
        --input=FILE/FOLDER    path of the padmet representing the network to convert
        --output=FILE/FOLDER    path of json output file
        --focus=FILE/STR    path of tabulated compounds/reaction/pathway or compounds name
        -v
"""
import docopt
import os
import json
import sys

from cobra.io import read_sbml_model
from collections import OrderedDict
from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def metexploreviz_export_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    input_file_folder = args["--input"]
    output_file = args["--output"]
    focus = args["--focus"]
    verbose = args["-v"]
    metexploreviz_export(input_file_folder, output_file, focus, verbose)


def create_json_from_padmet(input_file, verbose=False):
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
    nodes_in_json = {}
    #loading padmetSpec
    if verbose:
        print('Loading %s' %input_file)
    padmetSpec = PadmetSpec(input_file)

    nodes_data = {}
    links_data = {}
    for node in padmetSpec.dicOfNode.values():
        if node.type == 'reaction':
            if node.misc["DIRECTION"][0] == 'REVERSIBLE':
                reversibility = True
            else:
                reversibility = False
            pathways_ids = [rlt.id_out for rlt in padmetSpec.dicOfRelationIn[node.id] if rlt.type == "is_in_pathway"]
            if node.id not in nodes_in_json:
                nodes_in_json[node.id] = len(nodes_in_json)
                nodes_data[node.id] = {'name':node.id, 'id': nodes_in_json[node.id], 'reactionReversibility':reversibility, 'biologicalType':'reaction', 'selected': False, 'labelVisible': False, 'pathways': pathways_ids}
            for rlt in padmetSpec.dicOfRelationIn[node.id]:
                if rlt.type in ['consumes']:
                    reactant_id = rlt.id_out
                    reactant_compartment = rlt.misc['COMPARTMENT'][0]
                    if reactant_id not in nodes_in_json:
                        nodes_in_json[reactant_id] = len(nodes_in_json)
                        nodes_data[reactant_id] = {'name':reactant_id, 'id': nodes_in_json[reactant_id], 'biologicalType':'metabolite', 'compartment': reactant_compartment, 'selected': False, 'labelVisible': False, 'pathways': pathways_ids}
                    else:
                        nodes_data[reactant_id]['pathways'] = list(set(nodes_data[reactant_id]['pathways'] + pathways_ids))
                    json_id =  str(nodes_in_json[reactant_id]) + ' -- ' + str(nodes_in_json[node.id])
                    links_data[json_id] = {'id':json_id, 'source': nodes_in_json[reactant_id], 'target': nodes_in_json[node.id], 'interaction': 'in', 'reversible': str(reversibility)}

                if rlt.type in ['produces']:
                    product_id = rlt.id_out
                    product_compartment = rlt.misc['COMPARTMENT'][0]
                    if product_id not in nodes_in_json:
                        nodes_in_json[product_id] = len(nodes_in_json)
                        nodes_data[product_id] = {'name':product_id, 'id': nodes_in_json[product_id], 'biologicalType':'metabolite', 'biologicalType':'metabolite', 'compartment': product_compartment, 'selected': False, 'labelVisible': False, 'pathways': pathways_ids}
                    else:
                        nodes_data[product_id]['pathways'] = list(set(nodes_data[product_id]['pathways'] + pathways_ids))
                    json_id =  str(nodes_in_json[node.id]) + ' -- ' + str(nodes_in_json[product_id])
                    links_data[json_id] = {'id':json_id, 'source': nodes_in_json[product_id], 'target': nodes_in_json[node.id], 'interaction': 'in', 'reversible': str(reversibility)}

    return nodes_data, links_data


def create_json_from_sbml(input_file, verbose=False):
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
    nodes_in_json = {}
    #loading padmetSpec
    if verbose:
        print('Loading %s' %input_file)
    sbml_model = read_sbml_model(input_file)

    nodes_data = {}
    links_data = {}
    for rxn in sbml_model.reactions:
        rxn_id = rxn.id
        reversibility = rxn.reversibility

        if rxn_id not in nodes_in_json:
            nodes_in_json[rxn_id] = len(nodes_in_json)
            nodes_data[rxn_id] = {'name':rxn_id, 'id': nodes_in_json[rxn_id], 'reactionReversibility':reversibility, 'biologicalType':'reaction', 'selected': False, 'labelVisible': False}
        for reactant in rxn.reactants:
            reactant_id = reactant.id
            reactant_compartment = reactant.compartment
            if reactant_id not in nodes_in_json:
                nodes_in_json[reactant_id] = len(nodes_in_json)
                nodes_data[reactant_id] = {'name':reactant_id, 'id': nodes_in_json[reactant_id], 'biologicalType':'metabolite', 'compartment': reactant_compartment, 'selected': False, 'labelVisible': False}
            json_id =  str(nodes_in_json[reactant_id]) + ' -- ' + str(nodes_in_json[rxn_id])
            links_data[json_id] = {'id':json_id, 'source': nodes_in_json[reactant_id], 'target': nodes_in_json[rxn_id], 'interaction': 'in', 'reversible': str(reversibility)}

        for product in rxn.products:
            product_id = product.id
            product_compartment = product.compartment
            if product_id not in nodes_in_json:
                nodes_in_json[product_id] = len(nodes_in_json)
                nodes_data[product_id] = {'name':product_id, 'id': nodes_in_json[product_id], 'biologicalType':'metabolite', 'biologicalType':'metabolite', 'compartment': product_compartment, 'selected': False, 'labelVisible': False}
            json_id =  str(nodes_in_json[rxn_id]) + ' -- ' + str(nodes_in_json[product_id])
            links_data[json_id] = {'id':json_id, 'source': nodes_in_json[product_id], 'target': nodes_in_json[rxn_id], 'interaction': 'in', 'reversible': str(reversibility)}

    return nodes_data, links_data


def select_nodes_links(nodes_data, links_data, focus):
    keep_pathways = nodes_data[focus]['pathways']

    new_nodes_data = OrderedDict()
    for node_id in nodes_data:
        intersection_keep_pathways = list(set(nodes_data[node_id]['pathways']).intersection(set(keep_pathways)))

        if len(intersection_keep_pathways) > 0:
            new_nodes_data[node_id] = nodes_data[node_id]
            new_nodes_data[node_id]['pathways'] = intersection_keep_pathways

    new_node_ids = {str(new_nodes_data[node]['id']): str(index) for index, node in enumerate(new_nodes_data)}

    new_links_data = {}
    for link_id in links_data:
        node_id_1 = link_id.split(' -- ')[0]
        node_id_2 = link_id.split(' -- ')[1]
        if node_id_1 in new_node_ids and node_id_2 in new_node_ids:
            new_link_id = new_node_ids[node_id_1] + ' -- ' + new_node_ids[node_id_2]
            new_links_data[new_link_id] = links_data[link_id]
            new_links_data[new_link_id]['id'] = new_link_id
            new_links_data[new_link_id]['source'] = int(new_node_ids[node_id_1])
            new_links_data[new_link_id]['target'] = int(new_node_ids[node_id_2])

    return new_nodes_data, new_links_data


def metexploreviz_export(input_file_folder, output_file, focus, verbose=False):
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

    json_dicts = {}

    if input_type == 'file':
        _, file_extension = os.path.splitext(input_file_folder)
        json_dicts_tmp = {}
        json_dicts_tmp['nodes'] = []
        json_dicts_tmp['links'] = []
        if file_extension == '.padmet':
            nodes_data, links_data = create_json_from_padmet(input_file_folder, verbose)
        elif file_extension == '.sbml':
            nodes_data, links_data = create_json_from_sbml(input_file_folder, verbose)
        else:
            print('Incorrect format for {} Need .padmet of .sbml file.'.format(input_file_folder))
            sys.exit()
        if focus:
            nodes_data, links_data = select_nodes_links(nodes_data, links_data, focus)
        for node_id in nodes_data:
            json_dicts_tmp['nodes'].append(nodes_data[node_id])
        for link_id in links_data:
            json_dicts_tmp['links'].append(links_data[link_id])
        json_dicts.update(json_dicts_tmp)

    elif input_type == 'dir':
        for input_file in os.listdir(input_file_folder):
            json_dicts_tmp = {}
            json_dicts_tmp['nodes'] = []
            json_dicts_tmp['links'] = []
            input_path = os.path.join(input_file_folder, input_file)
            _, file_extension = os.path.splitext(input_path)
            if file_extension == '.padmet':
                nodes_data, links_data = create_json_from_padmet(input_path, verbose)
            elif file_extension == '.sbml':
                nodes_data, links_data = create_json_from_sbml(input_path, verbose)
            else:
                print('Incorrect format for {} Need .padmet of .sbml file.'.format(input_path))
                sys.exit()
            for node_id in nodes_data:
                json_dicts_tmp['nodes'].append(nodes_data[node_id])
            for link_id in links_data:
                json_dicts_tmp['links'].append(links_data[link_id])
            json_dicts.update(json_dicts_tmp)

    with open(output_file, 'w') as dumpfile:
        json.dump(json_dicts, dumpfile, indent=4)
