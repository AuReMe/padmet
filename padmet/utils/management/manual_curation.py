# -*- coding: utf-8 -*-
"""
Description:
    Update a padmetSpec by filling specific forms.

    1./ Create new reaction(s) to padmet file. \n
    - Get the template form with --template_new_rxn
    - Fill the template
    - set --data as path to the filled template

    2./ Add reaction(s) from padmetRef or remove reactions(s). \n
    - Get the template form with --template_add_delete_rxn
    - Fill the template
    - set --date as path to the filled template

    Update padmetSpec and create a new padmet (new_padmet) or overwrite the input

::

    usage:
        padmet manual_curation --padmetSpec=FILE --data=FILE [--padmetRef=FILE] [--output=FILE] [--tool=STR] [--category=STR] [-v]
        padmet manual_curation --template_new_rxn=FILE
        padmet manual_curation --template_add_delete_rxn=FILE

    option:
        -h --help    Show help.
        --padmetSpec=FILE    path to the padmet to update
        --padmetRef=FILE    path of the padmet representing the reference database
        --data=FILE    path to the form with data for curation
        --output=FILE    path to the output. if None. Overwriting padmetSpec
        --tool=STR    specification of the tool used to allow this curation: ex a tool of gapfilling (meneco)
        --category=STR    specification of the category of curation: ex if a reaction is added based on annotation info, use 'annotation'
        --template_new_rxn=FILE    create a form used to create new reaction, use this form as input for 'data' option
        --template_add_delete_rxn=FILE    create a form used to add or delete reaction, use this form as input for 'data' option
        -v    print info
"""
import csv
import docopt
import os

from padmet.classes import Relation, PadmetRef, PadmetSpec
from padmet.utils.sbmlPlugin import parseGeneAssoc


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def manual_curation_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    data_file = args["--data"]
    output = args["--output"]
    verbose = args["-v"]

    if data_file:
        if not os.path.exists(data_file):
            raise FileNotFoundError("No form curation file (--data/data_file) accessible at " + data_file)

        filename = os.path.splitext(os.path.basename(data_file))[0]
        source = filename

    category = args["--category"]
    tool = args["--tool"]
    if args["--template_new_rxn"]:
        output = args["--template_new_rxn"]
        template_new_rxn(output)
    elif args["--template_add_delete_rxn"]:
        output = args["--template_add_delete_rxn"]
        template_add_delete(output)
    else:
        padmetSpec = PadmetSpec(args["--padmetSpec"])
        if not output:
            output = args["--padmetSpec"]
        if args["--padmetRef"]:
            padmetRef = PadmetRef(args["--padmetRef"])
        else:
            padmetRef = None
        to_do = sniff_datafile(data_file)

        if to_do == "rxn_creator":
            rxn_creator(data_file, padmetSpec, output, padmetRef, source, tool, category, verbose)
        elif to_do == "add_delete_rxn":
            add_delete_rxn(data_file, padmetSpec, output, padmetRef, source, tool, category, verbose)


def sniff_datafile(data_file):
    """
    Read data_file and check which kind of data input it is.
    A reaction_creator file contains only 2 columns.
    Add reaction_add_delete more than 2. Basic, need to be improved.

    Parameters
    ----------
    data_file: str
        path to file of reaction_creator or reaction_add_delete.

    Returns
    -------
    str:
        "rxn_creator" or "add_delete_rxn"
    """
    with open(data_file, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        header = next(reader)
        if len(header) == 2:
            to_do = "rxn_creator"
        elif len(header) > 2:
            to_do = "add_delete_rxn"
        else:
            raise TypeError("Unable to read the file")
    return to_do

def rxn_creator(data_file, padmetSpec, output, padmetRef=None, source=None, tool=None, category="MANUAL", verbose=False):
    """
    Read a data_file (form created with template_new_rxn and filed), for each reaction
    to create, add the reaction in padmetSpec (only if the id of the reaction is not already in padmetSpec or in padmetRef if given)
    the source ensure the traceability of the reaction, its a simple tag ex 'pathway_XX_update'
    if not given the filename of data_file will be used.
    if a tool was used to infer the reaction, define tool='name_of_the_tool'
    the Padmet of reference padmetRef can be used to check that the reaction id is not
    already in the database and copy information from the database for existing compounds
    strongly recommended to give a padmetRef.

    Parameters
    ----------
    data_file: str
        path to file based on template_new_rxn()
    padmetSpec: padmet.classes.PadmetSpec
        padmet to update
    output: str
        path to the new padmet file
    source: str
        tag associated to the new reactions to create and add, used for traceability
    tool: str
        The eventual tool used to infer the reactions to create and add
    category: str
        The default category of the reaction added manually is 'MANUAL'. Must not be changed.
    padmetRef: padmet.classes.PadmetRef
        padmet containing the database of reference
    verbose: bool
        if True print information
    """
    if not source:
        filename = os.path.splitext(os.path.basename(data_file))[0]
        source = filename
    source = source.upper()
    if tool:
        tool = tool.upper()
    if not category:
        category = "MANUAL"


    dict_data = {}
    with open(data_file, 'r') as f:
        all_read = f.read()
    sep = csv.Sniffer().sniff(all_read).delimiter
    data = (line for line in all_read.splitlines() if len(line) != 0 and not line.startswith("#"))
    for line in data:
        #if len of value is 0 then TypeError raised
        try:
            attrib, value = line.split(sep)
        except TypeError:
            continue
        attrib = attrib.replace(" ", "")
        if attrib == "reaction_id":
            current_id = value
            dict_data[current_id] = {}
        else:
            try:
                dict_data[current_id][attrib] .append(value)
            except KeyError:
                dict_data[current_id][attrib] = [value]
    if verbose:
        print("%s reactions to add" %len(list(dict_data.keys())))
    for reaction_id, reaction_data in dict_data.items():
        if verbose:
            print("check if the id %s is already used" %reaction_id)
        if reaction_id in list(padmetSpec.dicOfNode.keys()):
            print("the id : %s is already associated to an other reaction in padmetSpec, choose an other" %reaction_id)
            continue
        if padmetRef is not None and reaction_id in list(padmetRef.dicOfNode.keys()):
            print("the id : %s is already associated to an other reaction in padmetRef, choose an other" %reaction_id)
            continue

        if verbose:
            print("Adding reaction %s" %reaction_id)
        reaction_rev = reaction_data["reversible"][0].lower()
        if reaction_rev.upper() == "TRUE":
            reaction_rev = "REVERSIBLE"
        elif reaction_rev.upper() == "FALSE":
            reaction_rev = "LEFT-TO-RIGHT"
        else:
            print("Please choose a value in ['true','false'] for the reversibility of the reaction: %s" %reaction_id)
            continue
        comment = reaction_data["comment"]
        node_misc = {"DIRECTION":[reaction_rev]}
        padmetSpec.createNode("reaction", reaction_id, node_misc)

        #reconstructionData:
        if tool:
            reconstructionData_id = reaction_id+"_reconstructionData_"+tool
            reconstructionData = {"SOURCE": [source], "CATEGORY":[category], "TOOL":[tool], "COMMENT":comment}
            if reconstructionData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                print("Warning: The reaction %s seems to be already added from the same source %s" %(reaction_id, tool))
        else:
            reconstructionData_id = reaction_id+"_reconstructionData_MANUAL"
            reconstructionData = {"SOURCE": [source], "CATEGORY":["MANUAL"], "COMMENT":comment}
            if reconstructionData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %reaction_id)

        reconstructionData_rlt = Relation(reaction_id, "has_reconstructionData", reconstructionData_id)
        padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])

        genes_assoc = reaction_data["linked_gene"][0]
        if genes_assoc:
            #suppData:
            if tool:
                suppData_id = reaction_id+"_SuppData_"+tool
                if suppData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(reaction_id, tool))
            else:
                suppData_id = reaction_id+"_SuppData_MANUAL"
                if suppData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %reaction_id)
            suppData = {"GENE_ASSOCIATION":[genes_assoc]}
            #create the node suppData and the relation has_suppData
            suppData_rlt = Relation(reaction_id, "has_suppData", suppData_id)
            padmetSpec.createNode("suppData", suppData_id, suppData, [suppData_rlt])

            all_genes = parseGeneAssoc(genes_assoc)
            nbGenes = len(all_genes)
            if verbose:
                print("%s is linked to %s genes" %(reaction_id, nbGenes))
            for gene_id in all_genes:
                try:
                    #check if gene already in the padmet
                    padmetSpec.dicOfNode[gene_id]
                except KeyError:
                    padmetSpec.createNode("gene", gene_id)
                #check if rxn already linked to gene x
                try:
                    linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[reaction_id]
                    if rlt.type == "is_linked_to" and rlt.id_out == gene_id][0]
                    #rxn already linked to gene x, update misc
                    try:
                        linked_rlt.misc["SOURCE:ASSIGNMENT"].append(source)
                    except KeyError:
                        linked_rlt.misc["SOURCE:ASSIGNMENT"] = [source]
                #rxn not linked to gene x
                except IndexError:
                    linked_rlt = Relation(reaction_id, "is_linked_to", gene_id, {"SOURCE:ASSIGNMENT":[source]})
                padmetSpec._addRelation(linked_rlt)

        if verbose:
            print("check if all metabolites are already in the network")
        try:
            for reactant_data in reaction_data["reactant"]:
                stoechio, metabo_id, compart = reactant_data.split(":")
                stoechio = stoechio.replace(",", ".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose:
                        print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose:
                                print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        padmetSpec.createNode("compound", metabo_id)
                        if verbose:
                            print(("new compound created: id = %s" %metabo_id))
                rlt = Relation(reaction_id, "consumes", metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio], "COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose: print("No reactants defined")

        try:
            for product_data in reaction_data["product"]:
                stoechio, metabo_id, compart = product_data.split(":")
                stoechio = stoechio.replace(",", ".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose:
                        print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose:
                                print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        padmetSpec.createNode("compound", metabo_id)
                        print("new compound created: id = %s" % metabo_id)
                rlt = Relation(reaction_id, "produces", metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio], "COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose:
                print("No products defined")
        if "pathway" in reaction_data.keys():
            pathways = reaction_data["pathway"][0].split(";")
            for pwy_id in pathways:
                try:
                    padmetSpec.dicOfNode[pwy_id]
                except KeyError:
                    if verbose:
                        print("%s not in the network" %pwy_id)
                    if padmetRef is not None:
                        if verbose:
                            print("Check if new pathway %s is in dbref" %pwy_id)
                        if pwy_id in padmetRef.dicOfNode.keys():
                            print("Warning the new pathway %s exist in the dbref, risk of overwritting data, change pwy id" %pwy_id)
                            continue
                    padmetSpec.createNode("pathway", pwy_id)
                    if verbose:
                        print(("new pathway created: id = %s" %pwy_id))
                rlt = Relation(reaction_id, "is_in_pathway", pwy_id)
                padmetSpec._addRelation(rlt)
    if verbose:
        print("Creating output: %s" % output)
    padmetSpec.generateFile(output)

def add_delete_rxn(data_file, padmetSpec, output, padmetRef=None, source=None, tool=None, category="MANUAL", verbose=False):
    """
    Read a data_file (form created with template_add_delete and filed), for each reaction
    if column 'Action' == 'add':
        add the reaction from padmetRef to padmetSpec.
    elif column 'Action' == 'delete':
        remove the reaction
    Can't add a reaction without a padmetRef !

    the source ensure the traceability of the reaction, its a simple tag ex 'pathway_XX_update'
    if not given the filename of data_file will be used.
    if a tool was used to infer the reaction, define tool='name_of_the_tool'

    Parameters
    ----------
    data_file: str
        path to file based on template_new_rxn()
    padmetSpec: padmet.classes.PadmetSpec
        padmet to update
    padmetRef: padmet.classes.PadmetRef
        padmet containing the database of reference
    output: str
        path to the new padmet file
    source: str
        tag associated to the new reactions to create and add, used for traceability
    tool: str
        The eventual tool used to infer the reactions to create and add
    category: str
        The default category of the reaction added manually is 'MANUAL'. Must not be changed.
    verbose: bool
        if True print information
    """
    if not source:
        filename = os.path.splitext(os.path.basename(data_file))[0]
        source = filename
    source = source.upper()
    if tool:
        tool = tool.upper()
    if not category:
        category = "MANUAL"

    with open(data_file, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        file_name = os.path.basename(data_file)
        file_name = os.path.splitext(file_name)[0]

        reader = csv.DictReader(csvfile, delimiter=dialect.delimiter)
        for row in reader:
            element_id, comment, action, genes_assoc = row["idRef"], row["Comment"], row["Action"], row.get("Genes", None)
            if action.upper() == "ADD":
                if padmetRef is None:
                    if verbose:
                        print("No given padmetRef, unable to copy %s" %element_id)
                else:
                    if verbose:
                        print("Adding: %s" %(element_id))
                    padmetSpec.copyNode(padmetRef, element_id)

                    #reconstructionData:
                    if tool:
                        reconstructionData_id = element_id+"_reconstructionData_"+tool
                        reconstructionData = {"SOURCE": [source], "CATEGORY":[category], "TOOL":[tool], "COMMENT":[comment]}
                        if reconstructionData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                            print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id, tool))
                    else:
                        reconstructionData_id = element_id+"_reconstructionData_MANUAL"
                        reconstructionData = {"SOURCE": [source], "CATEGORY":["MANUAL"], "COMMENT":[comment]}
                        if reconstructionData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                            print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %element_id)

                    reconstructionData_rlt = Relation(element_id, "has_reconstructionData", reconstructionData_id)
                    padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])

                    if genes_assoc:
                        #suppData:
                        if tool:
                            suppData_id = element_id+"_SuppData_"+tool
                            if suppData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                                print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id, tool))
                        else:
                            suppData_id = element_id+"_SuppData_MANUAL"
                            if suppData_id in list(padmetSpec.dicOfNode.keys()) and verbose:
                                print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %element_id)
                        suppData = {"GENE_ASSOCIATION":[genes_assoc]}
                        #create the node suppData and the relation has_suppData
                        suppData_rlt = Relation(element_id, "has_suppData", suppData_id)
                        padmetSpec.createNode("suppData", suppData_id, suppData, [suppData_rlt])

                        all_genes = parseGeneAssoc(genes_assoc)
                        nbGenes = len(all_genes)
                        if verbose:
                            print("%s is linked to %s genes" %(element_id, nbGenes))
                        for gene_id in all_genes:
                            try:
                                #check if gene already in the padmet
                                padmetSpec.dicOfNode[gene_id]
                            except KeyError:
                                padmetSpec.createNode("gene", gene_id)
                            #check if rxn already linked to gene x
                            try:
                                linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[element_id]
                                if rlt.type == "is_linked_to" and rlt.id_out == gene_id][0]
                                #rxn already linked to gene x, update misc
                                try:
                                    linked_rlt.misc["SOURCE:ASSIGNMENT"].append(source)
                                except KeyError:
                                    linked_rlt.misc["SOURCE:ASSIGNMENT"] = [source]
                            #rxn not linked to gene x
                            except IndexError:
                                linked_rlt = Relation(element_id, "is_linked_to", gene_id, {"SOURCE:ASSIGNMENT":[source]})
                            padmetSpec._addRelation(linked_rlt)

            elif action.upper() == "DELETE":
                if verbose:
                    print("deleting: %s" %(element_id))
                padmetSpec.delNode(element_id)
            elif action == "":
                print("Nothing to do for: %s" %(element_id))
            else:
                print("Action: %s unknown for %s" %(action, element_id))
                print("action must be = 'add' or 'delete' or ''")
                exit()
        padmetSpec.generateFile(output)


def template_new_rxn(output):
    """
    Generate template file used as input of rxn_creator function

    Parameters
    ----------
    output: str
        path for the template new_rxn to create
    """
    with open(output, 'w') as f:
        line = "\t".join(["reaction_id", "my_rxn"])+"\n"
        f.write(line)
        line = "\t".join(["comment", "reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible", "false"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene", "(gene_a or gene_b) and gene_c"])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product", "#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant", "1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant", "2.0:compound_b:c"])+"\n"
        f.write(line)
        line = "\t".join(["product", "1.0:compound_c:c"])+"\n"
        f.write(line)
        f.write("\n")
        line = "\t".join(["reaction_id", "my_rxn_2"])+"\n"
        f.write(line)
        line = "\t".join(["comment", "reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible", "true"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene", ""])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product", "#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant", "1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant", "2.0:compound_d:c"])+"\n"
        f.write(line)
        line = "\t".join(["product", "1.0:compound_c:c"])+"\n"
        f.write(line)

def template_add_delete(output):
    """
    Generate template file used as input of add_delete_rxn function

    Parameters
    ----------
    output: str
        path for the template rxn_add_delete to create
    """
    with open(output, 'w') as csvfile:
        fieldnames = ["idRef", "Comment", "Action", "Genes"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow({"idRef": 'rxn_id_1', 'Comment': 'Reaction deleted for x reason', "Genes":"", "Action":"delete"})
        writer.writerow({"idRef": 'rxn_id_2', 'Comment': 'Reaction added for x reason', "Genes":"(gene1 and gene2)", "Action":"add"})
        writer.writerow({"idRef": 'rxn_id_3', 'Comment': 'Reaction added for x reason', "Genes":"", "Action":"add"})

if __name__ == "__main__":
    main()
