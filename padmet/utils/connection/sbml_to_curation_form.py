# -*- coding: utf-8 -*-
"""
Description:
    extract 1 reaction (if rxn_id) or a list of reactions (if rxn_file) 
    from a sbml file to the form used in aureme for curation.
    For example use this script to extract specific missing reaction of a model to
    a just created metabolic network.

::

    usage:
        padmet sbml_to_curation_form --sbml=FILE --output=FILE --rxn_id=ID [--comment=STR] [--extract-gene] [-v]
        padmet sbml_to_curation_form --sbml=FILE --output=FILE --rxn_file=FILE [--comment=STR] [--extract-gene] [-v]

    options:
        -h --help     Show help.
        --sbml=FILE    path of the sbml.
        --output=FILE    form containing the reaction extracted, form used for manual curation in aureme.
        --rxn_id=FILE    id of one reaction to extract
        --rxn_file=FILE    file of reactions ids to extract, 1 id by line.
        --extract-gene    If true, extract also genes associated to reactions.
        --comment=STR    comment associated to the reactions in the form. Used to track sources of curation in aureme [default: "N.A"].
        -v   print info
"""
import docopt
import libsbml
import os

from padmet.utils.sbmlPlugin import convert_from_coded_id, parseNotes


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def sbml_to_curation_form_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    sbml_file = args["--sbml"]
    if args["--rxn_id"]:
        rxn_list = [args["--rxn_id"]]
    if args["--rxn_file"]:
        with open(args["--rxn_file"], 'r') as f:
            rxn_list = f.read().splitlines()
    output = args["--output"]
    comment = args["--comment"]
    verbose = args["-v"]
    extract_gene = args["--extract-gene"]
    sbml_to_curation(sbml_file, rxn_list, output, extract_gene=extract_gene, comment=comment, verbose=verbose)


def sbml_to_curation(sbml_file, rxn_list, output, extract_gene=False, comment="N.A", verbose=False):
    """
    Read a sbml file, check if each reaction ids are in the sbml, if no, raise ValueError
    Then create the form. this form can then be used with manual_curation.py

    Parameters
    ----------
    sbml_file: str
        path to sbml file
    rxn_list: list
        list of reaction id, ids must be identic as in the sbml, carrefull to encoded ids.
    output: str
        path to the form to create
    extract_gene: bool
        if true extract genes association
    comment: str
        Comment why the reaction will be added in the network for traceability.
    verbose: bool
        if True print information
    
    """
    if not os.path.exists(sbml_file):
        raise FileNotFoundError("No SBML file (--sbml/sbml_file) accessible at " + sbml_file)

    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    #check if reactions id are in model.
    if verbose:
        print("Check if reaction(s) are in sbml file")
    for rxn_id in rxn_list:
        if rxn_id in [r.id for r in listOfReactions]:
            if verbose:
                print("reaction %s found" %rxn_id)
        else:
            raise ValueError("/!\ reaction %s not found" %rxn_id)
            
    #create form output
    with open(output, 'w') as f:
        for rxn_id in rxn_list:
            rxn_sbml = listOfReactions.getElementBySId(rxn_id)
            rxn_id_decoded = convert_from_coded_id(rxn_id)[0]
            if verbose:
                print("extracting reaction %s, decoded id as %s" %(rxn_id, rxn_id_decoded))
            line = ["reaction_id",rxn_id_decoded]
            line = "\t".join(line)+"\n"
            f.write(line)
            line = ["comment",comment]
            line = "\t".join(line)+"\n"
            f.write(line)
            if rxn_sbml.reversible:
                line = ["reversible","true"]
            else:
                line = ["reversible","false"]
            line = "\t".join(line)+"\n"
            f.write(line)
            #check if have gene assoc
            if extract_gene:
                try:
                    gene_assoc = parseNotes(rxn_sbml)["GENE_ASSOCIATION"][0]
                    line = ["linked_gene", gene_assoc]
                except KeyError:
                    line = ["linked_gene", ""]
            else:
                    line = ["linked_gene", ""]

            line = "\t".join(line)+"\n"
            f.write(line)
            line = ["#reactant/product","#stoichio:compound_id:compart"]
            line = "\t".join(line)+"\n"
            f.write(line)
            reactants = rxn_sbml.getListOfReactants()
            products = rxn_sbml.getListOfProducts()
            for reactant in reactants:
                stoich = str(abs(reactant.getStoichiometry()))
                reactant_id,x, compart = convert_from_coded_id(reactant.getSpecies())
                line = ":".join([stoich, reactant_id, compart])
                line = "reactant"+"\t"+line+"\n"
                f.write(line)
            for product in products:
                stoich = str(abs(product.getStoichiometry()))
                product_id,x, compart = convert_from_coded_id(product.getSpecies())
                line = ":".join([stoich, product_id, compart])
                line = "product"+"\t"+line+"\n"
                f.write(line)
            f.write("\n")

