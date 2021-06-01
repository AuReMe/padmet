# -*- coding: utf-8 -*-
"""
Description:
    convert a padmet representing a metabolic network into a json compatible with MetExplore.

    usage:
        padmet metexploreviz_export --input=FILE --output=DIR [-v]

    options:
        -h --help     Show help.
        --input=FILE/FOLDER    path of the padmet representing the network to convert
        --output=FILE/FOLDER    path of json output file
        -v
"""
import docopt
import os
import json
import sys

from cobra.io import read_sbml_model
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
    verbose = args["-v"]
    metexploreviz_export(input_file_folder, output_file, verbose)


def create_json_from_padmet(input_file, verbose=False):
    """
    #TODO
    """
    json_dicts = {}
    json_dicts['nodes'] = []
    json_dicts['links'] = []

    nodes_in_json = {}
    #loading padmetSpec
    if verbose:
        print('Loading %s' %input_file)
    padmetSpec = PadmetSpec(input_file)
    padmetSpec_name = os.path.splitext(os.path.basename(input_file))[0]

    for node in padmetSpec.dicOfNode.values():
        if node.type == 'reaction':
            if node.misc["DIRECTION"][0] == 'REVERSIBLE':
                reversibility = True
            else:
                reversibility = False
            if node.id not in nodes_in_json:
                nodes_in_json[node.id] = len(nodes_in_json)
                json_dicts['nodes'].append({'name':node.id, 'id': nodes_in_json[node.id], 'reactionReversibility':reversibility, 'biologicalType':'reaction', 'selected': False, 'labelVisible': False})
            for rlt in padmetSpec.dicOfRelationIn[node.id]:
                if rlt.type in ['consumes']:
                    reactant_id = rlt.id_out
                    reactant_compartment = rlt.misc['COMPARTMENT'][0]
                    if reactant_id not in nodes_in_json:
                        nodes_in_json[reactant_id] = len(nodes_in_json)
                        json_dicts['nodes'].append({'name':reactant_id, 'id': nodes_in_json[reactant_id], 'biologicalType':'metabolite', 'compartment': reactant_compartment, 'selected': False, 'labelVisible': False})
                    json_id =  str(nodes_in_json[reactant_id]) + ' -- ' + str(nodes_in_json[node.id])
                    json_dicts['links'].append({'id':json_id, 'source': nodes_in_json[reactant_id], 'target': nodes_in_json[node.id], 'interaction': 'in', 'reversible': str(reversibility)})

                if rlt.type in ['produces']:
                    product_id = rlt.id_out
                    product_compartment = rlt.misc['COMPARTMENT'][0]
                    if product_id not in nodes_in_json:
                        nodes_in_json[product_id] = len(nodes_in_json)
                        json_dicts['nodes'].append({'name':product_id, 'id': nodes_in_json[product_id], 'biologicalType':'metabolite', 'biologicalType':'metabolite', 'compartment': product_compartment, 'selected': False, 'labelVisible': False})
                    json_id =  str(nodes_in_json[node.id]) + ' -- ' + str(nodes_in_json[product_id])
                    json_dicts['links'].append({'id':json_id, 'source': nodes_in_json[node.id], 'target': nodes_in_json[product_id], 'interaction': 'out', 'reversible': str(reversibility)})

    return json_dicts


def create_json_from_sbml(input_file, verbose=False):
    """
    #TODO
    """
    json_dicts = {}
    json_dicts['nodes'] = []
    json_dicts['links'] = []

    nodes_in_json = {}
    #loading padmetSpec
    if verbose:
        print('Loading %s' %input_file)
    sbml_model = read_sbml_model(input_file)

    for rxn in sbml_model.reactions:
        rxn_id = rxn.id
        reversibility = rxn.reversibility

        if rxn_id not in nodes_in_json:
            nodes_in_json[rxn_id] = len(nodes_in_json)
            json_dicts['nodes'].append({'name': rxn_id, 'id': nodes_in_json[rxn_id], 'reactionReversibility':reversibility, 'biologicalType':'reaction', 'selected': False, 'labelVisible': False})
        for reactant in rxn.reactants:
            reactant_id = reactant.id
            reactant_compartment = reactant.compartment
            if reactant_id not in nodes_in_json:
                nodes_in_json[reactant_id] = len(nodes_in_json)
                json_dicts['nodes'].append({'name':reactant_id, 'id': nodes_in_json[reactant_id], 'biologicalType':'metabolite', 'compartment': reactant_compartment, 'selected': False, 'labelVisible': False})
            json_id =  str(nodes_in_json[reactant_id]) + ' -- ' + str(nodes_in_json[rxn_id])
            json_dicts['links'].append({'id':json_id, 'source': nodes_in_json[reactant_id], 'target': nodes_in_json[rxn_id], 'interaction': 'in', 'reversible': str(reversibility)})

        for product in rxn.products:
            product_id = product.id
            product_compartment = product.compartment
            if product_id not in nodes_in_json:
                nodes_in_json[product_id] = len(nodes_in_json)
                json_dicts['nodes'].append({'name':product_id, 'id': nodes_in_json[product_id], 'biologicalType':'metabolite', 'compartment': product_compartment, 'selected': False, 'labelVisible': False})
            json_id =  str(nodes_in_json[rxn_id]) + ' -- ' + str(nodes_in_json[product_id])
            json_dicts['links'].append({'id':json_id, 'source': nodes_in_json[rxn_id], 'target': nodes_in_json[product_id], 'interaction': 'out', 'reversible': str(reversibility)})

    return json_dicts


def metexploreviz_export(input_file_folder, output_file, verbose):
    if os.path.isdir(input_file_folder):
        input_type = "dir"
    elif os.path.isfile(input_file_folder):
        input_type = "file"
    else:
        raise TypeError("%s is not a dir or a file or is not accessible." %(input_file_folder))

    json_dicts = {}

    if input_type == 'file':
        _, file_extension = os.path.splitext(input_file_folder)

        if file_extension == '.padmet':
            json_dicts_tmp = create_json_from_padmet(input_file_folder, verbose)
        elif file_extension == '.sbml':
            json_dicts_tmp = create_json_from_sbml(input_file_folder, verbose)
        else:
            print('Incorrect format for {} Need .padmet of .sbml file.'.format(input_file_folder))
            sys.exit()
        json_dicts.update(json_dicts_tmp)

    elif input_type == 'dir':
        for input_file in os.listdir(input_file_folder):
            input_path = os.path.join(input_file_folder, input_file)
            _, file_extension = os.path.splitext(input_path)
            if file_extension == '.padmet':
                json_dicts_tmp = create_json_from_padmet(input_path, verbose)
            elif file_extension == '.sbml':
                json_dicts_tmp = create_json_from_sbml(input_path, verbose)
            else:
                print('Incorrect format for {} Need .padmet of .sbml file.'.format(input_path))
                sys.exit()
            json_dicts.update(json_dicts_tmp)

    with open(output_file, 'w') as dumpfile:
        json.dump(json_dicts, dumpfile, indent=4)
