# -*- coding: utf-8 -*-
"""
Description:
    convert a padmet representing a database (padmetRef) and/or a padmet representing a model (padmetSpec)
    to tsv files for askomics.

    1./ Folder creation
    given the output directory. Create this directory if required and create a folder
    padmetRef filename and/or padmetSpec filename
    
    2./

    2.1/ For padmetRef:

        2.1.a/ Nodes
            get all reactions nodes => extract data from misc with extract_nodes(rxn_nodes, "reaction", "../rxn.tsv")

            get all compounds nodes => extract data from misc with extract_nodes(cpd_nodes, "compounds", "../cpd.tsv")

            get all pathways nodes => extract data from misc with extract_nodes(pwy_nodes, "pathway", "../pwy.tsv")

            get all xrefs nodes => extract data from misc with extract_nodes(xref_nodes, "xref", "../xref.tsv")
    
        2.1.b/ Relations
            for each rxn in rxn nodes:

                get all rlt consumes/produces => create list of data with extract_rxn_cpd(rxn_cpd_rlt)
                    fieldnames = "rxn_cpd","concerns@reaction","consumes@compound","produces@compound","stoichiometry","compartment"

                get all rlt is_in_pathway => create list of data with extract_rxn_pwy(rxn_pwy_rlt)
                    fieldnames = "rxn_pwy","concerns@reaction","in_pwy@pathway"

                get all rlt has_xref => create list of data with extract_entity_xref(rxn_xref_rlt)

            for each cpd in cpd nodes:

                get all rlt has_xref => update previous list of data with extract_entity_xref(cpd_xref_rlt)
                    fieldnames = "entity_xref","concerns@reaction","concerns@compound","has_xref@xref"

::

    usage:
        padmet padmet_to_tsv --padmetSpec=FILE [--padmetRef=FILE] --output_dir=DIR [-v]
        padmet padmet_to_tsv --padmetRef=FILE [--padmetSpec=FILE] --output_dir=DIR [-v]

    options:
        -h --help     Show help.
        --padmetSpec=FILE    path of the padmet representing the network to convert
        --padmetRef=FILE    path of the padmet representing the database
        --output_dir=DIR
        -v
"""
import docopt
import os
import csv

from padmet.classes import PadmetSpec, PadmetRef
from itertools import zip_longest


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def padmet_to_tsv_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    padmetSpec_file = args["--padmetSpec"]
    padmetRef_file = args["--padmetRef"]
    output_dir = args["--output_dir"]
    verbose = args["-v"]
    padmet_to_tsv(padmetSpec_file, padmetRef_file, output_dir, verbose)


def padmet_to_tsv(padmetSpec_file, padmetRef_file, output_dir, verbose=False):
    """
    #TODO
    """
    global all_rxn_nodes, all_cpd_nodes, all_pwy_nodes
    
    #check if output_dir exist, else create it
    if not output_dir.endswith("/"):
        output_dir += "/"
    if not os.path.isdir(output_dir):
        if verbose: print("Creating folder %s" %output_dir)
        os.makedirs(output_dir)
    #loading padmetSpec
    if padmetSpec_file:
        if verbose: print("Loading %s" %padmetSpec_file)
        padmetSpec = PadmetSpec(padmetSpec_file)
        padmetSpec_name = os.path.splitext(os.path.basename(padmetSpec_file))[0]
        padmetSpec_folder = output_dir+padmetSpec_name+"/"
        if not os.path.isdir(padmetSpec_folder):
            if verbose: print(("Creating folder %s" %padmetSpec_folder))
            os.makedirs(padmetSpec_folder)

    #if padmetRef given, create folder for padmetRef
    if padmetRef_file:
        if verbose: print("Loading %s" %padmetRef_file)
        padmetRef = PadmetRef(padmetRef_file)
        padmetRef_name = os.path.splitext(os.path.basename(padmetRef_file))[0]
        padmetRef_folder = output_dir+padmetRef_name+"/"
        if not os.path.isdir(padmetRef_folder):
            if verbose: print("Creating folder %s" %padmetRef_folder)
            os.makedirs(padmetRef_folder)

    #NODES
    #Converting nodes data to tsv format
    #Except for names nodes: adding them as attribu
    if padmetRef_file:
        if verbose: print("Extracting nodes from %s" %padmetRef_name)
        with open(padmetRef_folder+"metabolic_network.tsv", 'w') as f:
            fieldnames = ["metabolic_network","name"]
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            writer.writerow([padmetRef_name, padmetRef_name])

        if verbose: print("\tExtracting reactions")
        all_rxn_nodes = [node for node in list(padmetRef.dicOfNode.values()) if node.type == "reaction"]
        if all_rxn_nodes: extract_nodes(padmetRef, all_rxn_nodes, "reaction", padmetRef_folder+"rxn.tsv", {"in@metabolic_network":[padmetRef_name]})
        if verbose: print("\t%s reactions" %len(all_rxn_nodes))

        if verbose: print("\tExtracting compounds")
        all_cpd_nodes = set([padmetRef.dicOfNode[rlt.id_out] for rlt in padmetRef.getAllRelation() if rlt.type in ["consumes","produces"]])
        if all_cpd_nodes: extract_nodes(padmetRef, all_cpd_nodes, "compound", padmetRef_folder+"cpd.tsv")
        if verbose: print("\t%s compounds" %len(all_cpd_nodes))
    
        if verbose: print("\tExtracting pathways")
        all_pwy_nodes = [node for node in list(padmetRef.dicOfNode.values()) if node.type == "pathway"]
        if all_pwy_nodes: extract_nodes(padmetRef, all_pwy_nodes, "pathway", padmetRef_folder+"pwy.tsv")
        if verbose: print("\t%s pathways" %len(all_pwy_nodes))

        if verbose: print("\tExtracting xrefs")
        all_xrefs_nodes = [node for node in list(padmetRef.dicOfNode.values()) if node.type == "xref"]
        if all_xrefs_nodes: extract_nodes(padmetRef, all_xrefs_nodes, "xref", padmetRef_folder+"xref.tsv")
        if verbose: print("\t%s xrefs" %len(all_xrefs_nodes))

        #RELATIONS
        #Converting relations data to tsv format
        if verbose: print("Extracting relations from %s" %padmetRef_name)
        rxn_cpd_data = []
        rxn_pwy_data = []
        entity_xref_data = []
        if verbose: print("\tExtracting relations reaction-[consumes/produces]-compound")
        if verbose: print("\tExtracting relations reaction-is_inclued_in-pathway")
        if verbose: print("\tExtracting relations reactions-has_xref-xref")

        for rxn_node in all_rxn_nodes:
            rxn_id = rxn_node.id
            #if verbose: print("Reaction %s" %rxn_id)

            #all consumes/produces relations
            cp_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]]
            rxn_cpd_data += extract_rxn_cpd(cp_rlt)
            #all is_in_pathway relations
            pwy_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"]
            if pwy_rlt: rxn_pwy_data += extract_rxn_pwy(pwy_rlt)
            #all has_xref relations
            rxn_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type == "has_xref"]
            if rxn_xref_rlt: entity_xref_data += extract_entity_xref(rxn_xref_rlt, padmetRef)

        if verbose: print("\tExtracting relations compound-has_xref-xref")
        for cpd_node in all_cpd_nodes:
            cpd_id = cpd_node.id
            try:
                cpd_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[cpd_id] if rlt.type == "has_xref"]
                if cpd_xref_rlt: entity_xref_data += extract_entity_xref(cpd_xref_rlt, padmetRef)
            except KeyError:
                pass

        if verbose: print("\tExtracting relations pwy-has_xref-xref")
        for pwy_node in all_pwy_nodes:
            pwy_id = pwy_node.id
            try:
                pwy_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[pwy_id] if rlt.type == "has_xref"]
                if pwy_xref_rlt: entity_xref_data += extract_entity_xref(pwy_xref_rlt, padmetRef)
            except KeyError:
                pass
        
        if rxn_cpd_data:
            if verbose: print("\t\tCreating rxn_cpd.tsv")
            rxn_cpd_file(rxn_cpd_data, padmetRef_folder+"rxn_cpd.tsv")
        if rxn_pwy_data:
            if verbose: print("\t\tCreating rxn_pwy.tsv")
            rxn_pwy_file(rxn_pwy_data, padmetRef_folder+"rxn_pwy.tsv")
        if entity_xref_data:
            if verbose: print("\t\tCreating entity_xref.tsv")
            entity_xref_file(entity_xref_data, padmetRef_folder+"entity_xref.tsv")

    else:
        if verbose: print("No given padmetRef")
        all_rxn_nodes, all_cpd_nodes, all_pwy_nodes, all_xref_nodes=[[]]*4

    if padmetSpec_file:
        if verbose: print("Extracting nodes from %s" %padmetSpec_name)
        with open(padmetSpec_folder+"metabolic_network.tsv", 'w') as f:
            fieldnames = ["metabolic_network","name"]
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            writer.writerow([padmetSpec_name,padmetSpec_name])

        if verbose: print("\tExtracting reactions")
        spec_rxn_nodes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "reaction"]
        if all_rxn_nodes: extract_nodes(padmetSpec, spec_rxn_nodes, "reaction", padmetSpec_folder+"rxn.tsv", {"in@metabolic_network":[padmetSpec_name]})
        if verbose: print("\t%s reactions" %len(spec_rxn_nodes))

        if verbose: print("\tExtracting compounds")
        spec_cpd_nodes = set([padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.getAllRelation() if rlt.type in ["consumes","produces"]])
        if all_cpd_nodes: extract_nodes(padmetSpec, spec_cpd_nodes, "compound", padmetSpec_folder+"cpd.tsv")
        if verbose: print("\t%s compounds" %len(spec_cpd_nodes))
    
        if verbose: print("\tExtracting pathways")
        spec_pwy_nodes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "pathway"]
        if all_pwy_nodes: extract_nodes(padmetSpec, spec_pwy_nodes, "pathway", padmetSpec_folder+"pwy.tsv")
        if verbose: print("\t%s pathways" %len(spec_pwy_nodes))

        if verbose: print("\tExtracting xrefs")
        spec_xrefs_nodes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "xref"]
        if spec_xrefs_nodes: extract_nodes(padmetSpec, spec_xrefs_nodes, "xref", padmetSpec_folder+"xref.tsv")
        if verbose: print("\t%s xrefs" %len(spec_xrefs_nodes))

        if verbose: print("\tExtracting all genes")
        spec_genes_nodes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "gene"]
        if spec_genes_nodes: extract_nodes(padmetSpec, spec_genes_nodes, "gene", padmetSpec_folder+"gene.tsv", opt_col = {"in@metabolic_network":[padmetSpec_name]})
        if verbose: print("\t%s genes" %len(spec_genes_nodes))

        if verbose: print("\tExtracting all reconstructionData")
        spec_recData_nodes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "reconstructionData"]
        if spec_genes_nodes: extract_nodes(padmetSpec, spec_recData_nodes, "reconstructionData", padmetSpec_folder+"reconstructionData.tsv")
        if verbose: print("\t%s reconstructionData" %len(spec_recData_nodes))

        if verbose: print("Extracting relations from %s" %padmetSpec_name)
        rxn_cpd_data = []
        rxn_pwy_data = []
        rxn_gene_data = []
        entity_xref_data = []
        rxn_rec_data = []

        """
        fieldnames = ["rxn_reconstruction_info","concers@reaction","has_metadata@reconstruction_information","concerns@metabolic_network"]
        with open(padmetSpec_folder+"rxn_sources.tsv", 'w') as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            for rxn_node in spec_all_rxn_nodes:
                for src in rxn_node.misc.get("SOURCE",[]):
                    line = [rxn_node.id, src, padmetSpec_name]
                    line.insert(0,"_".join(line))
                    writer.writerow(line)            
        """
        if verbose: print("\tExtracting relations reaction-[consumes/produces]-compound")
        if verbose: print("\tExtracting relations reaction-is_in_pathway-pathway")
        if verbose: print("\tExtracting relations reactions-has_xref-xref")
        if verbose: print("\tExtracting relations reactions-is_linked_to-gene")
        if verbose: print("\tExtracting relations reactions-has_metadata-reconstructionData")

        for rxn_node in spec_rxn_nodes:
            rxn_id = rxn_node.id
            #if verbose: print("Reaction %s" %rxn_id)
            #all consumes/produces relations
            rxn_cpd_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]]
            rxn_cpd_data += extract_rxn_cpd(rxn_cpd_rlt)
            #all is_in_pathway relations
            rxn_pwy_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"]
            rxn_pwy_data += extract_rxn_pwy(rxn_pwy_rlt)
            #all has_xref relations
            rxn_xref_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "has_xref"]
            entity_xref_data += extract_entity_xref(rxn_xref_rlt, padmetSpec)
            #all is_linked_to relations
            rxn_gene_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]
            rxn_gene_data += extract_rxn_gene(rxn_gene_rlt)
            #all reconstructionData
            rxn_rec_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"]
            rxn_rec_data += extract_rxn_rec(rxn_rec_rlt)

        if verbose: print("\tExtracting relations compound-has_xref-xref")
        for cpd_node in spec_cpd_nodes:
            cpd_id = cpd_node.id
            #if verbose: print("Compound %s" %cpd_id)
            try:
                cpd_xref_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[cpd_id] if rlt.type == "has_xref"]
                entity_xref_data += extract_entity_xref(cpd_xref_rlt, padmetSpec)
            except KeyError:
                pass

        if padmetRef_file:
            if verbose: print("\tExtracting pathways's completion rate and creating pwy_rate.tsv")
            pwy_rate(padmetRef, padmetSpec, padmetSpec_name, padmetSpec_folder+"pwy_rate.tsv")
        else:
            if verbose: print("No padmetRef given unable to calculate pathway ratio")

        if rxn_cpd_data:
            if verbose: print("\t\tCreating rxn_cpd.tsv")
            rxn_cpd_file(rxn_cpd_data, padmetSpec_folder+"rxn_cpd.tsv")
        if rxn_pwy_data:
            if verbose: print("\t\tCreating rxn_pwy.tsv")
            rxn_pwy_file(rxn_pwy_data, padmetSpec_folder+"rxn_pwy.tsv")
        if entity_xref_data:
            if verbose: print("\t\tCreating entity_xref.tsv")
            entity_xref_file(entity_xref_data, padmetSpec_folder+"entity_xref.tsv")
        if rxn_gene_data:
            if verbose: print("\t\tCreating rxn_gene.tsv")
            rxn_gene_file(rxn_gene_data, padmetSpec_folder+"rxn_gene.tsv")
        if rxn_rec_data:
            if verbose: print("\t\tCreating rxn_reconstructionData.tsv")
            rxn_rec_file(rxn_rec_data, padmetSpec_folder+"rxn_reconstructionData.tsv")



def extract_nodes(padmet, nodes, entity_id, output, opt_col = {}):
    """
    for n nodes in nodes. for each node.misc = {A:['x'],B:['y','z']}
    create a file with line = [node.id,A[0],B[0]],[node.id,"",B[1]]
    the order is defined in fieldnames.
    merge common name and synonyms in 'name'
    
    #TODO
    """
    #set of keys from each misc dict
    all_keys = set()
    #opt_col: if wanted to add extra data
    [all_keys.update(list(node.misc.keys())) for node in nodes]
    if opt_col:
        fieldnames = [entity_id] + list(opt_col.keys()) + list(all_keys)
    else:
        fieldnames = [entity_id] + list(all_keys)
    
    with open(output, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        for node in nodes:
            try:
                names = [padmet.dicOfNode[rlt.id_out].misc["LABEL"] for rlt in padmet.dicOfRelationIn[node.id] if rlt.type == "has_name"]
                names = names[0]
                if "COMMON-NAME" in list(node.misc.keys()):
                    [node.misc["COMMON-NAME"].append(name) for name in names if name not in node.misc["COMMON-NAME"]]
                else:
                    node.misc["COMMON-NAME"] = names
            except (KeyError, IndexError):
                names = [""]
            if opt_col:
                node.misc.update(opt_col)
            data = [node.misc.get(col,[""]) for col in fieldnames]
            data = zip_longest(*data,fillvalue="")
            for d in data:
                d = list(d)
                d[0] = node.id
                writer.writerow(d)

def extract_rxn_cpd(rxn_cpd_rlt):
    """
    for rlt in rxn_cpd_rlt, append in data: [rxn_id,cpd_id(consumed),'',stoich,compartment]
    and/or [rxn_id,'',cpd_id(produced),stoich,compartment]. The value in index 0
    is a merge of all data to create a unique relation id
    
    #TODO
    """
    data = []
    for rlt in rxn_cpd_rlt:
        rxn_id = rlt.id_in
        cpd_id = rlt.id_out
        stoich = rlt.misc["STOICHIOMETRY"][0]
        if "n" in stoich: stoich = "1"
        if rlt.type == "consumes":
            line = [rxn_id, cpd_id, "", stoich, rlt.misc["COMPARTMENT"][0]]
            line.insert(0, "_".join(line))
        elif rlt.type == "produces":
            line = [rxn_id, "", cpd_id, stoich, rlt.misc["COMPARTMENT"][0]]
            line.insert(0, "_".join(line))
        data.append(line)
    return data

def rxn_cpd_file(data, output):
    """
    from data obtained with extract_rxn_cpd(), create file rxn_cpd
    
    #TODO
    """
    fieldnames = ["rxn_cpd","concerns@reaction","consumes@compound","produces@compound","stoichiometry","compartment"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

def extract_rxn_pwy(rxn_pwy_rlt):
    """
    for rlt in rxn_pwy_rlt, append in data: [rxn_id,pwy_id]. The value in index 0
    is a merge of all data to create a unique relation id
    
    #TODO
    """
    data = []
    for rlt in rxn_pwy_rlt:
        rxn_id = rlt.id_in
        pwy_id = rlt.id_out
        line = [rxn_id, pwy_id]
        line.insert(0,"_".join(line))
        data.append(line)
    return data

def extract_pwy(padmet):
    """
    from padmet return a dict, k = pwy_id, v = set of rxn_id in pwy
    
    #TODO
    """
    #get all pathways in ref in dict: k=pwy_id, v = set(rxn_id in pwy)
    #get all pathways in in spec dict: k=pwy_id, v = set(rxn_id in pwy)
    pathways_rxn_dict = dict([(node.id, set([rlt.id_in for rlt in padmet.dicOfRelationOut.get(node.id,[]) 
    if rlt.type == "is_in_pathway" and padmet.dicOfNode[rlt.id_in].type == "reaction"])) for node in list(padmet.dicOfNode.values()) if node.type == "pathway"])
    #pop k,v if len(v) == 0 (if not v)
    [pathways_rxn_dict.pop(k) for k,v in list(pathways_rxn_dict.items()) if not v]
    return pathways_rxn_dict

def pwy_rate(padmetRef, padmetSpec, metabolic_network, output):
    """
    pwy rate in padmetSpec is calculated based on padmetRef
    
    #TODO
    """
    all_pathways_dict = extract_pwy(padmetRef)    
    network_pathways_dict = extract_pwy(padmetSpec)
    fieldnames = ["pathway_rate","concerns@metabolic_network","concerns@pathway","rate"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        for pwy_id, all_rxns_set in list(all_pathways_dict.items()):
            #pwy_id = "PWY-5497"
            #rxns_set = network_pathways_dict["PWY-5497"]
            all_rxns_set = all_pathways_dict[pwy_id]
            nb_all_rxns = len(all_rxns_set)
            rxns_set = network_pathways_dict.get(pwy_id, set())
            nb_in_network = len(rxns_set.intersection(all_rxns_set))
            #rate nb rxn in network / nb total rxn in pwy
            rate = round(float(nb_in_network)/float(nb_all_rxns),2)
            rate = str(rate).replace(",",".")
            pwy_rate_id = "_".join([pwy_id, metabolic_network, rate])
            writer.writerow([pwy_rate_id, metabolic_network, pwy_id, rate])

def rxn_pwy_file(data, output):
    """
    #TODO
    """
    fieldnames = ["rxn_pwy","concerns@reaction","is_included_in@pathway"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]


def extract_entity_xref(entity_xref_rlt, padmet):
    """
    #TODO
    """
    data = []
    for rlt in entity_xref_rlt:
        entity_id = rlt.id_in
        xref_id = rlt.id_out
        if padmet.dicOfNode[entity_id].type == "reaction":
            line = [xref_id,entity_id, "", ""]
            line.insert(0, "_".join(line))
        elif padmet.dicOfNode[entity_id].type == "pathway":
            line = [xref_id, "", entity_id, ""]
            line.insert(0, "_".join(line))
        else:
            line = [xref_id, "", "", entity_id]
            line.insert(0, "_".join(line))
        data.append(line)
    return data

def entity_xref_file(data, output):
    """
    #TODO
    """
    fieldnames = ["entity_xref", "has_metadata@xref", "concerns@reaction","concerns@pathway","concerns@compound"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

def extract_rxn_gene(rxn_gene_rlt):
    """
    #TODO
    """
    data = []
    for rlt in rxn_gene_rlt:
        rxn_id = rlt.id_in
        gene_id = rlt.id_out
        all_src_data = rlt.misc["SOURCE:ASSIGNMENT"]
        for src_data in all_src_data:
            try:
                src, assignment = src_data.split(":")
            except ValueError:
                src = src_data
                assignment = ""
            recData_id = rxn_id + "_reconstructionData_" + src
            line = [rxn_id, gene_id, recData_id]
            line.insert(0, "_".join(line))
            line.append(assignment)
            data.append(line)
    return data

def rxn_gene_file(data, output):
    """
    #TODO
    """
    fieldnames = ["rxn_gene_recData","concerns@reaction","is_linked_to@gene","was_obtained_with@reconstructionData", "assignment_evidence"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

def extract_rxn_rec(rxn_rec_rlt):
    """
    #TODO
    """
    data = []
    for rlt in rxn_rec_rlt:
        rxn_id = rlt.id_in
        recData_id = rlt.id_out
        line = [rxn_id, recData_id]
        data.append(line)
    return data

def rxn_rec_file(data, output):
    """
    #TODO
    """
    fieldnames = ["reaction","has_metadata@reconstructionData"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]


