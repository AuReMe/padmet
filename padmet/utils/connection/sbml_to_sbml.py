#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create sbml from sbml. Use it to change sbml level.

"""

import os
import sys

from datetime import datetime
from multiprocessing import Pool
from padmet.classes import PadmetSpec
from padmet.utils.connection import sbmlGenerator

def run_sbml_to_sbml(multiprocess_data):
    """Turn sbml to sbml.

    Parameters
    ----------
    multiprocess_data: dict
        pathname to species sbml file, pathname to output sbml file, new sbml level

    Returns
    -------
    bool:
        True if sbml file exists
    """
    padmet = sbml_to_padmet(sbml=multiprocess_data['sbml_file'], db=None, version=None, source_tool=None,
                            source_category=None, source_id=None, mapping=None, verbose=None)
    sbmlGenerator.padmet_to_sbml(padmet, multiprocess_data['sbml_output_file'], sbml_lvl=multiprocess_data['new_sbml_level'], verbose=False)

    if multiprocess_data['sbml_output_file'] and not os.access(multiprocess_data['sbml_output_file'], os.W_OK):
        try:
            open(multiprocess_data['sbml_output_file'], 'w').close()
            os.unlink(multiprocess_data['sbml_output_file'])
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True

def sbml_to_padmet(sbml, db, version, source_tool, source_category, source_id, mapping, verbose):
    db='NA'
    version='NA'
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")

    if os.path.isdir(sbml):
        sbml_type = "dir"
    elif os.path.isfile(sbml):
        sbml_type = "file"
    else:
        raise TypeError("%s is not a dir or a file" %(sbml))
    padmet_to_update = PadmetSpec()

    POLICY_IN_ARRAY = [['class','is_a_class','class'], ['class','has_name','name'], ['class','has_xref','xref'], ['class','has_suppData','suppData'],
                    ['compound','is_a_class','class'], ['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                    ['gene','is_a_class','class'], ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                    ['pathway','is_a_class','class'], ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'],
                    ['protein','is_a_class','class'], ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                    ['protein','is_in_species','class'],
                    ['reaction','is_a_class','class'], ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],
                    ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
    dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":db,"version":version}}
    padmet_to_update.setInfo(dbNotes)
    padmet_to_update.setPolicy(POLICY_IN_ARRAY)

    #if sbml is a directory, recover all file path in a list. if no => only one file: create a list with only this file
    sbml_mapping_dict = {}
    if sbml_type == "dir":
        path = sbml
        if not path.endswith("/"):
            path += "/"
        all_files = [i for i in next(os.walk(path))[2] if not i.startswith(".~lock")]
        for sbml_file in [i for i in all_files if i.endswith(".sbml") or i.endswith(".xml")]:
            mapping_file = os.path.splitext(sbml_file)[0] + "_dict.csv"
            if mapping_file not in all_files:
                mapping_file = None
            else:
                mapping_file = path+mapping_file
            sbml_file = path+sbml_file
            sbml_mapping_dict[sbml_file] = mapping_file
    else:
        sbml_mapping_dict[sbml] = mapping

    for sbml_file, mapping_file in list(sbml_mapping_dict.items()):
        if mapping_file:
            force = False
        else:
            force = True

        padmet_to_update.updateFromSbml(sbml_file, verbose = verbose, source_category = source_category, source_id = source_id,
                                        source_tool = source_tool, mapping_file = mapping_file, force = force )

    if len(list(sbml_mapping_dict.keys())) == 0:
        if verbose: print("No sbml found in %s" %sbml)

    return padmet_to_update

def from_sbml_to_sbml(input_sbml, output_sbml, new_sbml_level, cpu=1):
    """Turn sbml to sbml.

    Parameters
    ----------
    input_sbml: str
        pathname to species sbml file/folder
    output_sbml: str
        pathname to output sbml file/folder
    new_sbml_level: int
        new sbml level
    cpu: int
        number of cpu

    Returns
    -------
    str:
        pathname to output sbml file/folder
    """
    try:
        new_sbml_level = int(new_sbml_level)
    except:
        print('SBML level must be an int.')

    if new_sbml_level not in [2, 3]:
        sys.exit('New SBML level must be 2 or 3.')

    sbml_to_sbml_pool = Pool(processes=int(cpu))

    multiprocess_datas = []
    if os.path.isdir(input_sbml):
        if not os.path.isdir(output_sbml):
            try:
                os.makedirs(output_sbml)
            except OSError:
                print('Impossible to create output folder.')
        for sbml_file in os.listdir(input_sbml):
            multiprocess_data = {}
            multiprocess_data['sbml_file'] = input_sbml + '/' + sbml_file
            multiprocess_data['sbml_output_file'] = output_sbml + '/' + sbml_file
            multiprocess_data['new_sbml_level'] = new_sbml_level
            multiprocess_datas.append(multiprocess_data)

    elif os.path.isfile(input_sbml):
        multiprocess_data = {}
        multiprocess_data['sbml_file'] = input_sbml
        multiprocess_data['sbml_output_file'] = output_sbml
        multiprocess_data['new_sbml_level'] = new_sbml_level
        multiprocess_datas.append(multiprocess_data)

    sbml_checks = sbml_to_sbml_pool.map(run_sbml_to_sbml, multiprocess_datas)

    sbml_to_sbml_pool.close()
    sbml_to_sbml_pool.join()

    if all(sbml_checks):
        return output_sbml

    else:
        print("Error during sbml creation.")

