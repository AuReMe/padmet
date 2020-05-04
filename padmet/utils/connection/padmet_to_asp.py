#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Convert PADMet to ASP following these predicats:
    common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)
    direction(reaction_id, reaction_direction). reaction_direction in[LEFT-TO-RIGHT,REVERSIBLE]
    ec_number(reaction_id, ec(x,x,x)).
    catalysed_by(reaction_id, enzyme_id).
    uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
    in_pathway(reaction_id, pathway_id).
    reactant(reaction_id, compound_id, stoechio_value).
    product(reaction_id, compound_id, stoechio_value).
    is_a(compound_id, class_id).
    is_a(pathway_id, pathway_id).

::

    usage:
        padmet padmet_to_asp --padmet=FILE --output=FILE [-v]

    option:
        -h --help     Show help.
        --padmet=FILE    path to padmet file to convert.
        --output=FILE    path to output file in lp format.
        -v    print info.
"""
import docopt

from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def padmet_to_asp_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    padmet_file = args["--padmet"]
    output = args["--output"]
    verbose = args["-v"]

    padmet_to_asp(padmet_file, output, verbose)


def padmet_to_asp(padmet_file, output, verbose = False):
    """
    Convert PADMet to ASP following these predicats:
    common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)
    direction(reaction_id, reaction_direction). reaction_direction in[LEFT-TO-RIGHT,REVERSIBLE]
    ec_number(reaction_id, ec(x,x,x)).
    catalysed_by(reaction_id, enzyme_id).
    uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
    in_pathway(reaction_id, pathway_id).
    reactant(reaction_id, compound_id, stoechio_value).
    product(reaction_id, compound_id, stoechio_value).
    is_a(compound_id, class_id).
    is_a(pathway_id, pathway_id).

    Parameters
    ----------
    padmet_file: str
        the path to padmet file to convert
    output: str
        the path to the output to create
    verbose: bool
        print informations
    """

    padmet = PadmetSpec(padmet_file)

    with open(output, 'w') as f:
        #recover all reactions's data
        reactions = [node for node in padmet.dicOfNode.values()
        if node.type == "reaction"]
        nb_reactions = len(reactions)

        #set where to store compounds and enzymes already recovered
        allCompounds_enzymes = set()
        count = 0
        for reaction_node in reactions:
            count += 1
            reaction_id = reaction_node.id
            if verbose: print("Reaction "+str(count)+"/"+str(nb_reactions)+"\t"+reaction_id)
            
            #common_name
            try:
                reaction_name = reaction_node.misc["COMMON_NAME"][0]
                line = asp_synt("common_name",[reaction_id, reaction_name])
                line += "\n"
                f.write(line)
            except KeyError:
                pass
            #direction
            try:
                reaction_direction = reaction_node.misc["DIRECTION"][0]
            except KeyError:
                reaction_direction = "REVERSIBLE"
            line = asp_synt("direction",[reaction_id, reaction_direction])
            line += "\n"
            f.write(line)
            
            #EC-Number
            try:
                ECs = reaction_node.misc["EC_NUMBER"]
                for ec_number in ECs:
                    ec_number = ec_number.replace("EC-","").replace(".",",")
                    line = "ec_number"+"(\""+reaction_id+"\", ec("+ec_number+"))."
                    line += "\n"
                    f.write(line)
            except KeyError: pass
            
            #Enzymes
            try:
                enzymes = [rlt.id_in for rlt in padmet.dicOfRelationOut[reaction_id]            
                if rlt.type == "catalyses"]
                for enzyme_id in enzymes:
                    #check if enzyme has uniprot_id
                    if enzyme_id not in allCompounds_enzymes:
                        allCompounds_enzymes.add(enzyme_id)
                        try:
                            uniprot_ref = [padmet.dicOfNode[rlt.id_out].misc["ID"][0]
                            for rlt in padmet.dicOfRelationIn[enzyme_id]
                            if rlt.type == "has_xref"
                            and ("UNIPROT" in padmet.dicOfNode[rlt.id_out].misc["DB"][0])] 
                            
                            for uniprot_id in uniprot_ref:
                                line = asp_synt("uniprotID",[enzyme_id, uniprot_id])
                                line += "\n"
                                f.write(line)
                        except KeyError: pass
                        try:
                            enzyme_name = padmet.dicOfNode[enzyme_id].misc["COMMON_NAME"][0]
                            line = asp_synt("common_name",[enzyme_id, enzyme_name])
                            line += "\n"
                            f.write(line)
                        except KeyError:
                            pass
                    line = asp_synt("catalysed_by",[reaction_id, enzyme_id])
                    line += "\n"
                    f.write(line)
            except KeyError: pass
            
            #pathways
            pathways = [rlt.id_out for rlt in padmet.dicOfRelationIn[reaction_id]
            if rlt.type == "is_in_pathway"]
            for pathway_id in pathways:
                line = asp_synt("in_pathway",[reaction_id, pathway_id])
                line += "\n"
                f.write(line)
        
            #reactants
            reactants_data = [(rlt.id_out, rlt.misc["STOICHIOMETRY"][0]) 
            for rlt in padmet.dicOfRelationIn[reaction_id] if rlt.type == "consumes"]
            for reactant_id, stoechio in reactants_data:
                try:
                    int(stoechio)
                except ValueError:
                    stoechio = "1"
                line = asp_synt("reactant",[reaction_id, reactant_id])
                line = line[:-2]+","+stoechio+line[-2:]+"\n"
                f.write(line)
                #class of reactant*
                if reactant_id not in allCompounds_enzymes:
                    allCompounds_enzymes.add(reactant_id)
                    try:
                        reactant_classes = [rlt.id_out for rlt in padmet.dicOfRelationIn[reactant_id]
                        if rlt.type == "is_a_class"]
                        for class_id in reactant_classes:
                            line = asp_synt("is_a",[reactant_id, class_id])
                            line += "\n"
                            f.write(line)
                    except KeyError: pass
                    try:
                        reactant_name = padmet.dicOfNode[reactant_id].misc["COMMON_NAME"][0]
                        line = asp_synt("common_name",[reactant_id, reactant_name])
                        line += "\n"
                        f.write(line)
                    except KeyError:
                        pass
            #products
            products_data = [(rlt.id_out, rlt.misc["STOICHIOMETRY"][0]) 
            for rlt in padmet.dicOfRelationIn[reaction_id] if rlt.type == "produces"]
            for product_id, stoechio in products_data:
                try:
                    int(stoechio)
                except ValueError:
                    stoechio = "1"
                line = asp_synt("product",[reaction_id, product_id])
                line = line[:-2]+","+stoechio+line[-2:]+"\n"
                f.write(line)
                #class of product
                if product_id not in allCompounds_enzymes:
                    allCompounds_enzymes.add(product_id)
                    try:
                        product_classes = [rlt.id_out for rlt in padmet.dicOfRelationIn[product_id]
                        if rlt.type == "is_a_class"]
                        for class_id in product_classes:
                            line = asp_synt("is_a",[product_id, class_id])
                            line += "\n"
                            f.write(line)
                    except KeyError: pass
                    try:
                        product_name = padmet.dicOfNode[product_id].misc["COMMON_NAME"][0]
                        line = asp_synt("common_name",[product_id, product_name])
                        line += "\n"
                        f.write(line)
                    except KeyError:
                        pass
        #recover all pathway's data
        pathways = [node for node in padmet.dicOfNode.values()
        if node.type == "pathway"]
        nb_pathways = len(pathways)
        count = 0
        for pathway_node in pathways:
            count += 1
            pathway_id = pathway_node.id
            if verbose: print("Pathway "+str(count)+"/"+str(nb_pathways)+"\t"+pathway_id)
            
            #common_name
            try:
                pathway_name = padmet.dicOfNode[pathway_id].misc["COMMON_NAME"][0]
                line = asp_synt("common_name",[pathway_id, pathway_name])
                line += "\n"
                f.write(line)
            except KeyError:
                pass
            is_in_pathway = [rlt.id_out for rlt in padmet.dicOfRelationIn[pathway_id]
            if rlt.type == "is_in_pathway"]
            for sub_pathway_id in is_in_pathway:
                line = asp_synt("is_a",[pathway_id, sub_pathway_id])
                line += "\n"
                f.write(line)

def asp_synt(pred, list_args):
    """
    create a predicat for asp
    
    example: asp_synt("direction",["R1","REVERSIBLE"]) => "direction('R1','reversible')."

    Parameters
    ----------
    pred: str
        the predicat
    list_args: list
        list of atoms to put in the predicat

    Returns
    -------
    str
        the predicat 'pred(''list_args[0]'',''list_args[1]'',...,''list_args[n]'').'
    """
    list_args = ["\""+str(i)+"\"" for i in list_args]
    rez = str(pred)+"("+",".join(list_args)+")."
    return rez


