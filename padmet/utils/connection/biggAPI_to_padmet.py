# -*- coding: utf-8 -*-
"""
Description:
    Require internet access !

    Allows to extract the bigg database from the API to create a padmet.

    1./ Get all reactions universal id from http://bigg.ucsd.edu/api/v2/universal/reactions, escape reactions of biomass.

    2./ Using async_list, extract all the informations for each reactions (compounds, stochio, name ...)

    3./ Need to use sleep time to avoid to lose the server access.

    4./ Because the direction fo the reaction is not set by default in bigg. 
    We get all the models where the reaction is and the final direction will the one found
    in more than 75%

    5./ Also extract xrefs

::

    usage:
        padmet biggAPI_to_padmet --output=FILE [--pwy_file=FILE] [-v]

    options:
        -h --help     Show help.
        --output=FILE    path to output, the padmet file.
        --pwy_file=FILE   add kegg pathways from pathways file, line:'pwy_id, pwy_name, x, rxn_id'.
        -v   print info.
"""
try:
    from gevent import monkey as curious_george
except ImportError:
    raise ImportError('Requires gevent, requests and grequests, try:\npip install gevent requests grequests')

curious_george.patch_all(thread=False, select=False)

from padmet.classes import Relation, instantiate_padmet

import docopt
import os

try:
    import requests
except ImportError:
    raise ImportError('Requires gevent, requests and grequests, try:\npip install gevent requests grequests')

try:
    import grequests
except ImportError:
    raise ImportError('Requires gevent, requests and grequests, try:\npip install gevent requests grequests')


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def biggAPI_to_padmet_cli(command_args):
    #parsing args
    args = docopt.docopt(__doc__, argv=command_args)
    output = args["--output"]
    verbose = args["-v"]
    pwy_file = args["--pwy_file"]
    biggAPI_to_padmet(output, pwy_file, verbose)


def biggAPI_to_padmet(output, pwy_file=None, verbose=False):
    """
    Extract BIGG database using the api. Create a padmet file.
    Escape reactions of biomass.
    Require internet access !

    Allows to extract the bigg database from the API to create a padmet.

    1./ Get all reactions universal id from http://bigg.ucsd.edu/api/v2/universal/reactions, escape reactions of biomass.
    2./ Using async_list, extract all the informations for each reactions (compounds, stochio, name ...)
    3./ Need to use sleep time to avoid to lose the server access.
    4./ Because the direction fo the reaction is not set by default in bigg. 
    We get all the models where the reaction is and the final direction will the one found
    in more than 75%
    5./ Also extract xrefs

    Parameters
    ----------
    output: str
        path to output, the padmet file.
    pwy_file: str
        path to pathway file, add kegg pathways, line:'pwy_id, pwy_name, x, rxn_id'.
    verbose: bool
        if True print information
    """
    padmet_id = os.path.splitext(os.path.basename(output))[0]
    padmetRef = instantiate_padmet("PadmetRef", None, padmet_id, "BIGG", "1.5", verbose)

    list_of_relation = []
    if verbose: print("Getting all reactions ids")
    url_bigg = 'http://bigg.ucsd.edu/api/v2/'
    raw_data = requests.get(url_bigg + "universal/reactions").json()['results']
    all_reactions_ids = [rxn_dict['bigg_id'] for rxn_dict in raw_data if not rxn_dict['bigg_id'].startswith("BIOMASS")]
    if verbose: print("%s reactions to extract" %(len(all_reactions_ids)))

    """
    if verbose: print("Extracting informations... Wait")
    step = 100
    rxn_lower_index = -(step)
    rxn_upper_index = 0
    rxn_responses = []
    all_range = len(all_reactions_ids)/step

    for i in range(all_range):
        async_list = []
        rxn_lower_index += step
        rxn_upper_index += step

        for rxn_id in all_reactions_ids[rxn_lower_index:rxn_upper_index]:
            action_item = grequests.get(url_bigg + "universal/reactions/" +rxn_id)
            async_list.append(action_item) 
        new_responses = [r.json() for r in grequests.map(async_list)]
        rxn_responses += new_responses
        print("%s/%s done" %(len(rxn_responses),len(all_reactions_ids)))

    if rxn_upper_index != len(all_reactions_ids):
        async_list = []
        last_index = len(all_reactions_ids) - rxn_upper_index
        rxn_lower_index += step
        rxn_upper_index += last_index
        for rxn_id in all_reactions_ids[rxn_lower_index:rxn_upper_index]:
            action_item = grequests.get(url_bigg + "universal/reactions/" +rxn_id)
            async_list.append(action_item) 
        new_responses = [r.json() for r in grequests.map(async_list)]
        rxn_responses += new_responses
    """
    if verbose: print("updating padmet")
    count = 0
    all_reactions_ids = [i for i in all_reactions_ids if 'biomass' not in i.upper()]
    for rxn_id in [i for i in all_reactions_ids if not i.startswith("BIOMASS")]:
        count += 1
        if verbose: print("reaction: %s, %s/%s" %(rxn_id, count, len(all_reactions_ids)))
        if rxn_id not in list(padmetRef.dicOfNode.keys()):
            rxn_response = requests.get(url_bigg + "universal/reactions/" +rxn_id)
            rxn_dict = rxn_response.json()
    
    
            rxn_metabolites = rxn_dict["metabolites"]
            if len(rxn_metabolites) > 1:
                rxn_id = rxn_dict['bigg_id']
                rxn_name = rxn_dict["name"]
            
                all_models_id = [i["bigg_id"] for i in rxn_dict["models_containing_reaction"]]
                async_list = []
                for model_id in all_models_id:
                    action_item = grequests.get(url_bigg + "models/"+ model_id +"/reactions/"+ rxn_id)
                    async_list.append(action_item)  
                models_responses = [r.json() for r in grequests.map(async_list)]
                all_lower_bound = [i["results"][0]["lower_bound"] for i in models_responses]
                ratio_not_rev = float(all_lower_bound.count(0))/float(len(all_lower_bound))
                if verbose: print("Reaction not reversible in %s/%s model(s)" %(all_lower_bound.count(0), len(all_lower_bound)))
                if ratio_not_rev >= 0.75:
                    rxn_direction = "LEFT-TO-RIGHT"
                    if verbose: print("Reaction not reversible")
                else:
                    rxn_direction = "REVERSIBLE"
                    if verbose: print("Reaction reversible")
                padmetRef.createNode("reaction",rxn_id,{"COMMON_NAME":[rxn_name],"DIRECTION":[rxn_direction]})
            
                rxn_xrefs = rxn_dict["database_links"]
    
                xref_id = rxn_id+"_xrefs"
                xref_node = padmetRef.createNode("xref", xref_id)
                has_xref_rlt = Relation(rxn_id, "has_xref", xref_id)
                list_of_relation.append(has_xref_rlt)
    
                for db, k in list(rxn_xrefs.items()):
                    _id = k[0]["id"]
                    if db in list(xref_node.misc.keys()) and _id not in xref_node.misc[db]:
                        xref_node.misc[db].append(_id)
                    else:
                        xref_node.misc[db] = [_id]
    
                for metabo_dict in rxn_metabolites:
                    metabo_id = metabo_dict["bigg_id"]
                    metabo_name = metabo_dict["name"]
                    metabo_compart = metabo_dict["compartment_bigg_id"]
                    metabo_stoich = metabo_dict["stoichiometry"]
                    try:
                        padmetRef.dicOfNode[metabo_id]
                    except KeyError:
                        padmetRef.createNode("compound",metabo_id,{"COMMON_NAME":[metabo_name]})
                    if metabo_stoich < 0:
                        consumes_rlt = Relation(rxn_id,"consumes",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                        list_of_relation.append(consumes_rlt)
                    else:
                        produces_rlt = Relation(rxn_id,"produces",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                        list_of_relation.append(produces_rlt)
        else: 
            if verbose: print("%s already in padmet" %rxn_id)
            continue                                        
    if verbose: print("Adding all relations")
    count = 0
    for rlt in list_of_relation:
        count += 1
        if verbose: print("relation %s/%s" %(count, len(list_of_relation)))
        try:
            padmetRef.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmetRef.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmetRef.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmetRef.dicOfRelationOut[rlt.id_out] = [rlt]
    
    if pwy_file:
        if not os.path.exists(pwy_file):
            raise FileNotFoundError("No KEGG Pathway file (--pwy_file/pwy_file) accessible at " + pwy_file)

        add_kegg_pwy(pwy_file, padmetRef, verbose)

    if verbose: print("Generating file: %s" %output)
    padmetRef.generateFile(output)

def add_kegg_pwy(pwy_file, padmetRef, verbose = False):
    """
    #TODO
    """
    global list_of_relation
    with open(pwy_file, 'r') as f:
        for data in [line.split("\t") for line in f.read().splitlines()][1:]:
            pwy_id, name, ec, rxn_id = data
            try:
                pwy_node = padmetRef.dicOfNode[pwy_id]
            except KeyError:
                pwy_node = padmetRef.createNode("pathway", pwy_id)
            if name:
                try:
                    pwy_node.misc["COMMON_NAME"].append(name)
                except KeyError:
                    pwy_node.misc["COMMON_NAME"] = [name]
            if rxn_id:
                if rxn_id in list(padmetRef.dicOfNode.keys()):
                    pwy_rlt = Relation(rxn_id,"is_in_pathway",pwy_id)
                    padmetRef._addRelation(pwy_rlt)
                else:
                    if verbose: print("%s in pwy %s but not in padmet" %(rxn_id, pwy_id))
            
