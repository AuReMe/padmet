# -*- coding: utf-8 -*-
"""
Description:
    The standard output of meneco return ids of reactions corresponding to the solution for gapfilling.

    The ids are those from the sbml and so they are encoded.

    This script extract the solution corresponding to the union of reactions
    "Computing union of reactions from all completion"
    Based on padmetRef return a file with more information for each reaction.

    ex: RXN__45__5

    RXN-5, common_name, ec-number, Formula (with id),Formula (with cname),Action,Comment
    Also, the output can be used as input of the script update_padmetSpec.py
    In the column Action: 'add' => To add the reaction, '' => to do nothing

    Comment: the reason of adding the reaction (ex: added for gap-filling by meneco)

::

    usage:
        padmet enhanced_meneco_output --meneco_output=FILE --padmetRef=FILE --output=FILE [-v]

    options:
        -h --help     Show help.
        --meneco_output=FILE    pathname of a meneco run' result
        --padmetRef=FILE    path to padmet file corresponding to the database of reference (the repair network)
        --output=FILE    path to tsv output file
"""
import docopt
import os

from padmet.classes import PadmetRef
from padmet.utils.sbmlPlugin import get_all_decoded_version


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def enhanced_meneco_output_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    meneco_output_file = args["--meneco_output"]
    output = args["--output"]
    verbose = args["-v"]
    padmetRef = PadmetRef(args["--padmetRef"])
    enhanced_meneco_output(meneco_output_file, padmetRef, output, verbose)


def enhanced_meneco_output(meneco_output_file, padmetRef, output, verbose=False):
    """
    The standard output of meneco return ids of reactions corresponding to the solution for gapfilling.
    The ids are those from the sbml and so they are encoded.
    This script extract the solution corresponding to the union of reactions
    "Computing union of reactions from all completion"
    Based on padmetRef return a file with more information for each reaction.

    ex: RXN__45__5
    RXN-5, common_name, ec-number, Formula (with id),Formula (with cname),Action,Comment
    Also, the output can be used as input for manual_curation
    In the column Action: 'add' => To add the reaction, '' => to do nothing
    Comment: the reason of adding the reaction (ex: added for gap-filling by meneco)
    
    Parameters
    ----------
    meneco_output_file: str
        pathname of a meneco run' result
    padmetRef: padmet.padmetRef
        path to padmet file corresponding to the database of reference (the repair network)
    output: str
        path to tsv output file
    verbose: bool
        if True print information    
    """
    if not os.path.exists(meneco_output_file):
        raise FileNotFoundError("No Meneco result (--meneco_output/meneco_output_file) accessible at " + meneco_output_file)

    with open(meneco_output_file,'r') as f:
        #recovering union reactions
        file_in_array = f.read().splitlines()
        start_index = None
        for line in file_in_array:
            if line.startswith("Computing union of reactions from all completion"):
                start_index = file_in_array.index(line) + 1
        #recover reactions, delete ' " ' and space.
        if start_index is None:
            print("No line starting with: Computing union of reactions from all completion. Enable to extracts reactions")
            #return
        encoded_reactions = [line.strip().replace("\"","") 
        for line in file_in_array[start_index:]]
        nb_reactions = len(encoded_reactions)
    if verbose: print("%s reactions to check" %nb_reactions)
    with open(output,'w') as f:
        header = ["idRef","Common name","EC-number","Formula (with id)","Formula (with cname)","Action","Comment", "Genes"]
        header = "\t".join(header)+"\n"
        f.write(header)
        for encoded_id in encoded_reactions:
            decoded_reactions = get_all_decoded_version(encoded_id, "reaction")
            reaction_id = set(decoded_reactions).intersection(set(padmetRef.dicOfNode.keys()))
            if reaction_id:
                reaction_id = list(reaction_id)[0]
                reac_node = padmetRef.dicOfNode[reaction_id]
                try:
                    ec = reac_node.misc["EC_NUMBER"][0]
                except KeyError:
                    ec = "Unknown"
                try:
                    common_name = reac_node.misc["COMMON_NAME"][0]
                except KeyError:
                    common_name = "Unknown"

                direction = reac_node.misc["DIRECTION"][0]
                if direction == "REVERSIBLE":
                    direction = " <=> "
                elif direction == "LEFT-TO-RIGHT":
                    direction = " => "
                else:
                    direction = " =>/<=> "

                id_reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out+"["+rlt.misc["COMPARTMENT"][0]+"]"
                for rlt in padmetRef.dicOfRelationIn.get(reaction_id,None) 
                if rlt.type == "consumes"]
                id_products = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out+"["+rlt.misc["COMPARTMENT"][0]+"]"
                for rlt in padmetRef.dicOfRelationIn.get(reaction_id,None) 
                if rlt.type == "produces"]
                idRef_formula = " + ".join(id_reactants)+direction+" + ".join(id_products)
                
                try:
                    cname_reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetRef.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0]+"["+rlt.misc["COMPARTMENT"][0]+"]" 
                    for rlt in padmetRef.dicOfRelationIn.get(reaction_id,None) 
                    if rlt.type == "consumes"]
                    cname_products = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetRef.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0]+"["+rlt.misc["COMPARTMENT"][0]+"]"
                    for rlt in padmetRef.dicOfRelationIn.get(reaction_id,None)
                    if rlt.type == "produces"]
                    cname_formula = " + ".join(cname_reactants)+direction+" + ".join(cname_products)
                except KeyError:
                    cname_formula = ""
                    
                line = [reaction_id, common_name, ec, idRef_formula, cname_formula, "add", "Added for gapfilling", ""]
                line = "\t".join(line)+"\n"
                f.write(line)
            else:
                print("%s not found in padmetRef" %reaction_id)
                pass

