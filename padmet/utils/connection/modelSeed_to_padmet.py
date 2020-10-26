# -*- coding: utf-8 -*-
"""
Description:
    From modelSeed reactions and pathways files creates a padmet file.

::

    usage:
        padmet modelSeed_to_padmet --output=FILE --rxn_file=FILE --pwy_file=FILE [-v]

    options:
        -h --help     Show help.
        --output=FILE    path of the padmet file to create
        --rxn_file=FILE   path to json file of modelSeed reactions
        --pwy_file=FILE   path to pathway reactions association from modelSeed
        -v   print info.
"""
import docopt
import csv
import json
import os

from padmet.classes import Relation, instantiate_padmet


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def modelSeed_to_padmet_cli(command_args):
    #parsing args
    args = docopt.docopt(__doc__, argv=command_args)
    output = args["--output"]
    verbose = args["-v"]
    rxn_file = args["--rxn_file"]
    pwy_file = args["--pwy_file"]
    modelSeed_to_padmet(rxn_file, pwy_file, output, verbose)


def modelSeed_to_padmet(rxn_file, pwy_file, output, verbose=False):
    """
    #TODO
    """
    global list_of_relation
    padmet_id = os.path.splitext(os.path.basename(output))[0]
    padmetRef = instantiate_padmet("PadmetRef", None, padmet_id, "MODELSEED", "1.0", verbose)

    list_of_relation = []

    if not os.path.exists(rxn_file):
        raise FileNotFoundError("No json file of modelSeed reactions (--rxn_file/rxn_file) accessible at " + rxn_file)

    if not os.path.exists(pwy_file):
        raise FileNotFoundError("No pathway reactions association file from modelSeed (--pwy_file/pwy_file) accessible at " + pwy_file)

    rxn_data = json.load(open(rxn_file))
    #remove biomass rxn:
    rxn_data.pop("rxn12985")
    if verbose: print("updating padmet")
    count = 0
    for rxn_id, rxn_dict in list(rxn_data.items()):
        count += 1
        if verbose: print("reaction: %s, %s/%s" %(rxn_id, count, len(rxn_data)))
        try:
            if not rxn_dict["compound_ids"]:
                raise KeyError
        except KeyError:
            print(rxn_id)
            continue
        if rxn_id not in list(padmetRef.dicOfNode.keys()):
            if rxn_dict["reversibility"] == ">":
                rxn_direction = "LEFT-TO-RIGHT"
            else:
                rxn_direction = "REVERSIBLE"
            rxn_name = rxn_dict["name"]
            padmetRef.createNode("reaction",rxn_id,{"COMMON_NAME":[rxn_name],"DIRECTION":[rxn_direction]})
                            
            rxn_metabolites = rxn_dict["stoichiometry"].split(";")
    
            for metabo_data in rxn_metabolites:
                metabo_data = metabo_data.replace("???","\"")
                try:
                    metabo_temp, metabo_name = metabo_data.split("\"")[:2]
                    metabo_stoich, metabo_id, metabo_compart = metabo_temp.split(":")[:3]
                except ValueError:
                    metabo_stoich, metabo_id, metabo_compart, metabo_name = metabo_data.split(":")[:4]
                    
                metabo_stoich = float(metabo_stoich)
                #from modelSeed github
                if metabo_compart == "0":
                    metabo_compart = "c"
                elif metabo_compart == "1":
                    metabo_compart = "e"
                elif metabo_compart == "2":
                    metabo_compart = "p"
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
    with open(pwy_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t' )
        pwy_raw_data = [row for row in reader]
    for pwy_raw in pwy_raw_data:
        pwy_id = pwy_raw["Source ID"]
        pwy_names = [pwy_raw["Name"],pwy_raw["Aliases"]]
        rxn_ids = pwy_raw["Reactions"].split("|")
        try:
            padmetRef.dicOfNode[pwy_id]
        except KeyError:
            padmetRef.createNode("pathway",pwy_id,{"COMMON_NAME":pwy_names})
        for rxn_id in rxn_ids:
            pwy_rlt = Relation(rxn_id,"is_in_pathway",pwy_id)
            list_of_relation.append(pwy_rlt)
        
        
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
    
    """
    if pwy_file:
        add_kegg_pwy(pwy_file, padmetRef, verbose)
    """
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
    padmetRef.generateFile("/home/maite/Documents/data/bigg/bigg_v2.padmet")
            
    

