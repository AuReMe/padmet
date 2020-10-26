#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Contains all necessary functions to generate wikiPages from a padmet file and update 
    a wiki online. Require WikiManager module (with wikiMate,Vendor)

::

    usage:
        padmet wikiGenerator --padmet=FILE/DIR --output=DIR --wiki_id=STR [--database=STR] [--padmetRef=FILE] [--log_file=FILE] [-v]
        padmet wikiGenerator --aureme_run=DIR --padmetSpec=ID -v

    options:
        -h --help     Show help.
        --padmet=FILE    path to padmet file.
        --output=DIR    path to folder to create with all wikipages in subdir.
        --wiki_id=STR    id of the wiki.
        --padmetRef=FILE    path to padmet of reference, ex: metacyc_xx.padmet, if given, able to calcul pathway rate completion.
        --log_file=FILE    log file from an aureme run, use this file to create a wikipage with all the command used during the aureme run.
        --aureme_run=DIR    can use an aureme run as input, will use from config file information for model_id and log_file and padmetRef.
        -v    print info.
"""
import docopt
import os
import shutil
import re

from padmet.classes import PadmetRef, PadmetSpec
from itertools import chain
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def wikiGenerator_cli(command_args):
    #files to upload: folder genomic_data, all sbml in output ortho, annot, external, seeds, targets
    args = docopt.docopt(__doc__, argv=command_args)
    padmet = args["--padmet"]
    verbose = args["-v"]
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    wiki_id = args["--wiki_id"]
    output = args["--output"]
    log_file = args["--log_file"]
    database = args["--database"]

    wikiGenerator(padmet, output, wiki_id, padmetRef, database, log_file, verbose)


def wikiGenerator(padmet, output, wiki_id, padmetRef=None, database=None, log_file=None, verbose=False):
    """
    #TODO
    """
    #if padmet is a directory: padmetFiles = list of padmet path, a padmet file must end with .padmet extension
    #padmetFiles = [/path/padmet_a.padmet, /path/padmet_b.padmet]
    #else: padmet could be many path separated by ';', ex: padmet  = "/path_a/padmet_a.padmet;/path_b/padmet_b.padmet"
    #padmetFiles = [/path_a/padmet_a.padmet, /path_b/padmet_b.padmet]
    if os.path.isdir(padmet):
        path = padmet
        all_files = [i for i in next(os.walk(path))[2] if not i.startswith(".~lock")]
        padmetFiles = [os.path.join(path, i) for i in all_files if i.endswith(".padmet")]
        if len(padmetFiles) == 0:
            raise IOError("No padmet found in %s" %path)
        
    else:
        padmetFiles = padmet.split(";")
    
    #if padmetRef given, extract all pathways associated to reactions (dont take in account pathways linked to pathway)
    #create a dict global_pwy_rxn_dict with k = pathway id, v = set of rxn_id
    global_pwy_rxn_dict = dict()
    if padmetRef:
        for rlt in [rlt for rlt in padmetRef.getAllRelation() if rlt.type == "is_in_pathway" and padmetRef.dicOfNode[rlt.id_in].type == "reaction"]:
            pwy_id = rlt.id_out
            rxn_id = rlt.id_in
            try:
                global_pwy_rxn_dict[pwy_id].add(rxn_id)
            except KeyError:
                global_pwy_rxn_dict[pwy_id] = set([rxn_id])

    #if database given, will add external links to specific database.
    #For example if database == metacyc, then will add external link to page title to metacyc website
    #Work only for 'metacyc' or 'bigg'
    if database: database = database.lower()
    if database == "metacyc":
        ext_link = {"reaction": "http://metacyc.org/META/NEW-IMAGE?object=",
                    "pathway": "http://metacyc.org/META/NEW-IMAGE?object=",
                    "metabolite": "http://metacyc.org/META/NEW-IMAGE?object="}
    elif database == "bigg":
        ext_link = {"reaction":"http://bigg.ucsd.edu/universal/reactions/",
                    "metabolite":"http://bigg.ucsd.edu/universal/metabolites/",
                    "pathway":"http://www.genome.jp/dbget-bin/www_bget?"}
    else:
        ext_link = {}

    #create all folder archi in output folder
    createDirectory(output, verbose)
    #total_padmet_data will store all the information from one to n padmet to finally use this data to create the pages
    """
    total_padmet_data["reaction"] = {'rxn-1':{"basic_attrib": dict(),
                                             "formula":dict(),
                                             "padmet_source": set(),
                                             "gene_assoc": dict(),
                                             "padmet_gene": dict(),
                                             "pathway_assoc":dict(),
                                             "reconstruction_data":dict(),
                                             "xref":dict(),
                                             "misc": dict()}}
    total_padmet_data["gene"] = {'gene-1':{"basic_attrib": dict(),
                                           "padmet_source":set(),
                                           "reac_assoc": dict(),
                                           "pathway_assoc":dict()}
    total_padmet_data["pathway"] = {'pwy-1':{"basic_attrib": dict(),
                                             "padmet_source":set(),
                                             "reac_assoc": dict(),
                                             "completion_rate":dict(),
                                             "total_reac_assoc":set()}
    total_padmet_data["metabolite"] = {'cpd-1':{"basic_attrib": dict(),
                                                "consumed_by": dict(),
                                                "produced_by": dict(),
                                                "unknown_dir_by": dict()}
    total_padmet_data["padmet_source"] = {"reaction":dict(),
                                          "gene": dict(),
                                          "pathway": dict(),
                                          "misc":dict()}
    total_padmet_data["reconstruction"] = {"category":dict(),
                                           "tool": dict(),
                                           "comment": dict(),
                                           "source": dict()}
    """
    total_padmet_data = {"reaction":{},
                         "gene":{},
                         "pathway":{},
                         "metabolite":{},
                         "padmet_source":{},
                         "reconstruction":{"category":dict(), "tool": dict(), "comment": dict(), "source": dict()},
                         "misc":{}
                         }
    #For each padmet file path in padmetFiles, extract all data with extract_padmet_data and add them to the dict total_padmet_data
    for padmetFile in padmetFiles:
        if verbose: print("Reading %s" %padmetFile)
        total_padmet_data = extract_padmet_data(padmetFile, total_padmet_data, global_pwy_rxn_dict, padmetRef, verbose)

    #Use reduce_padmet_data function to compress information. For example if all reactions have same values, remplace the list of padmet to tags
    total_padmet_data = reduce_padmet_data(total_padmet_data, verbose)

    #Create each categories pages
    for rxn_id, rxn_dict_data in total_padmet_data["reaction"].items():
        output_file = os.path.join(output, "reactions", rxn_id.replace("/", "__47__"))
        create_biological_page("reaction", rxn_id, rxn_dict_data, total_padmet_data, ext_link, output_file, padmetRef, verbose)
    for gene_id, gene_dict_data in total_padmet_data["gene"].items():
        output_file = os.path.join(output, "genes", gene_id)
        create_biological_page("gene", gene_id, gene_dict_data, total_padmet_data, ext_link, output_file, padmetRef, verbose)
    for pwy_id, pwy_dict_data in total_padmet_data["pathway"].items():
        output_file = os.path.join(output, "pathways", pwy_id.replace("/", "__47__"))
        create_biological_page("pathway", pwy_id, pwy_dict_data, total_padmet_data, ext_link, output_file, padmetRef, verbose)
    for cpd_id, cpd_dict_data in total_padmet_data["metabolite"].items():
        output_file = os.path.join(output, "metabolites", cpd_id)
        create_biological_page("metabolite", cpd_id, cpd_dict_data, total_padmet_data, ext_link, output_file, padmetRef, verbose)
    for org_id, org_dict_data in total_padmet_data["padmet_source"].items():
        output_file = os.path.join(output, "organisms", org_id)
        create_biological_page("organism", org_id, org_dict_data, total_padmet_data, ext_link, output_file, padmetRef, verbose)
    #Create navigation pages
    navigation_folder = os.path.join(output, "navigation")
    create_navigation_page(total_padmet_data, navigation_folder, verbose)

    if log_file:
        create_log_page(log_file, os.path.join(output, "navigation"))

    if len(padmetFiles) == 1:
        venn_path = os.path.join(output,"files","venn.png")
        create_venn(total_padmet_data, venn_path, verbose)
        main_page_path = os.path.join(output,"navigation","Main_Page")
        create_main(total_padmet_data, wiki_id, main_page_path, verbose)

def extract_padmet_data(padmetFile, total_padmet_data, global_pwy_rxn_dict=None, padmetRef=None, verbose=False):
    """
    total_padmet_data: k in ['reaction', 'gene', 'organism', 'pathway', ...]
    if k = 'reaction', v = {'misc':{},'gene_assoc':}
    
    For reaction in padmetFile:
        if reaction_id not in total_padmet_data["reaction"].keys():
            add total_padmet_data["reaction"][reaction_id][padmet_source] = dict()
            
    else, add data only if differents from first
    
    #TODO
    """
    #name of current padmet
    padmet_source = os.path.basename(padmetFile).replace(".padmet","")
    total_padmet_data["padmet_source"][padmet_source] = {"reaction":dict(), "gene": dict(), "pathway": dict(), "misc":dict()}
    current_padmet_dict = total_padmet_data["padmet_source"][padmet_source]
    padmet = PadmetSpec(padmetFile)

    all_rxns = [node for node in padmet.dicOfNode.values() if node.type == "reaction"]
    count = 0
    
    for rxn_node in all_rxns:
        count += 1
        rxn_id = rxn_node.id
        #rxn_id = '1-ACYLGLYCEROL-3-P-ACYLTRANSFER-RXN'
        #rxn_node = padmet.dicOfNode[rxn_id]
        if verbose:
            print("Extracting reaction:%s/%s\t%s" %(count, len(all_rxns), rxn_id), end='\r')
            
        current_padmet_dict["reaction"][rxn_id] = dict()

        #if first occurrenc, dict instantiation, in total_padmet_data["reaction"]
        if rxn_id not in total_padmet_data["reaction"].keys():
            total_padmet_data["reaction"][rxn_id] = {"basic_attrib": dict(), "formula":dict(), "padmet_source": set(), "gene_assoc": dict(), "padmet_gene": dict(), "pathway_assoc":dict(), "reconstruction_data":dict(), "xref":dict(), "misc": dict()}
            #If padmetRef given, search for external links ex ID in bigg, kegg, etc...
            if padmetRef:
                try:
                    xref_node = [padmetRef.dicOfNode[rlt.id_out] for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type == "has_xref"][0]
                    total_padmet_data["reaction"][rxn_id]["xref"] = dict(xref_node.misc)
                except (IndexError, KeyError) as e:
                    pass

        current_rxn_dict = total_padmet_data["reaction"][rxn_id]
        current_rxn_dict["padmet_source"].add(padmet_source)
        current_rxn_dict["reconstruction_data"][padmet_source] = list()

        #Use update_basic_attrib to extract data from node.misc and add them to 'basic_attrib'
        current_rxn_dict = update_basic_attrib(rxn_node, current_rxn_dict, padmet_source)                    
        list_names = [padmet.dicOfNode[rlt.id_out].misc["LABEL"] for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "has_name"]
        if list_names:
            current_names = set()
            [current_names.update(i) for i in list_names]
            if "synonymous" not in current_rxn_dict["basic_attrib"].keys():
                current_rxn_dict["basic_attrib"]["synonymous"] = dict()
            for name in current_names:
                name = name.lower()
                try:
                    current_rxn_dict["basic_attrib"]["synonymous"][name].add(padmet_source)
                except KeyError:
                    current_rxn_dict["basic_attrib"]["synonymous"][name] = set([padmet_source])

        #get all relations where reaction consumes or produces X
        consume_produce_rlts = [rlt for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]]
        for rlt in consume_produce_rlts:
            cpd_id = rlt.id_out
            cpd_node = padmet.dicOfNode[cpd_id]
            #For each compound, update total_padmet_data["metabolite"]
            if cpd_id not in total_padmet_data["metabolite"].keys():
                total_padmet_data["metabolite"][cpd_id] = {"basic_attrib": dict(), "consumed_by": dict(), "produced_by": dict(), "unknown_dir_by": dict()}
            #extract basic attrib of cpd
            current_cpd_dict = total_padmet_data["metabolite"][cpd_id]
            current_cpd_dict = update_basic_attrib(cpd_node, current_cpd_dict, padmet_source)                    

            #store compound depending on the direction of reaction
            if rxn_node.misc['DIRECTION'][0] == "LEFT-TO-RIGHT":
                direction_str = " '''=>''' "
                if rlt.type == "consumes":
                    try:
                        current_cpd_dict['consumed_by'][rxn_id].add(padmet_source)
                    except KeyError:
                        current_cpd_dict['consumed_by'][rxn_id] = set([padmet_source])
                elif rlt.type == "produces":
                    try:
                        current_cpd_dict['produced_by'][rxn_id].add(padmet_source)
                    except KeyError:
                        current_cpd_dict['produced_by'][rxn_id] = set([padmet_source])

            elif rxn_node.misc['DIRECTION'][0] == "REVERSIBLE":
                direction_str = " '''<=>''' "
                for by_type in ["consumed_by", "produced_by"]:
                    try:
                        current_cpd_dict[by_type][rxn_id].add(padmet_source)
                    except KeyError:
                        current_cpd_dict[by_type][rxn_id] = set([padmet_source])
            else:
                direction_str = " '''=>/<=>''' "
                try:
                    current_cpd_dict['unknown_dir_by'][rxn_id].add(padmet_source)
                except KeyError:
                    current_cpd_dict['unknown_dir_by'][rxn_id] = set([padmet_source])
            

        # Recovering the formula
        #reactants and products: list of str: ['stoichio [[cpd_id]][cpd_compart]',...]
        reactants = set([(rlt.misc["STOICHIOMETRY"][0], rlt.id_out, rlt.misc["COMPARTMENT"][0])
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "consumes"])
        products = set([(rlt.misc["STOICHIOMETRY"][0], rlt.id_out, rlt.misc["COMPARTMENT"][0])
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "produces"])

        reactants_id_str = [reac_tuple[0]+" [["+reac_tuple[1]+"]]["+reac_tuple[2]+"]"
                            for reac_tuple in  sorted(reactants, key=lambda reac_tuple: reac_tuple[1])]
        products_id_str = [prod_tuple[0]+" [["+prod_tuple[1]+"]]["+prod_tuple[2]+"]"
                            for prod_tuple in  sorted(products, key=lambda prod_tuple: prod_tuple[1])]
        #join each list with '+' and direction in center
        formula_id = " '''+''' ".join(reactants_id_str)+direction_str+" '''+''' ".join(products_id_str)
        try:
            current_rxn_dict["formula"][formula_id].add(padmet_source)
        except KeyError:
            current_rxn_dict["formula"][formula_id] = set([padmet_source])
            
        #get all pathways assocaited to the reaction    
        pathways_ids = [rlt.id_out for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"]
        for pwy_id in pathways_ids:
            pwy_node = padmet.dicOfNode[pwy_id]
            current_padmet_dict["pathway"][pwy_id] = dict()

            #for each pathway update total_padmet_data["pathway"]
            if pwy_id not in current_rxn_dict["pathway_assoc"].keys():
                current_rxn_dict["pathway_assoc"][pwy_id] = dict()
            if pwy_id not in total_padmet_data["pathway"].keys():
                total_padmet_data["pathway"][pwy_id] = {"basic_attrib": dict(), "padmet_source":set(), "reac_assoc": dict(), "completion_rate":dict(), "total_reac_assoc":set()}
            current_pwy_dict = total_padmet_data["pathway"][pwy_id]
            current_pwy_dict = update_basic_attrib(pwy_node, current_pwy_dict, padmet_source)      
            current_pwy_dict["padmet_source"].add(padmet_source)
            if rxn_id not in current_pwy_dict["reac_assoc"].keys():
                current_pwy_dict["reac_assoc"][rxn_id] = set()

            #recovering the nb of reactions associated to the pathway and the completion rate if global_pwy_rxn_dict given (if padmetRef given)
            total_reac_assoc = set()
            nb_total_reac_assoc = "n.a"
            if global_pwy_rxn_dict:
                if pwy_id in global_pwy_rxn_dict.keys():
                    total_reac_assoc = global_pwy_rxn_dict[pwy_id]
                    nb_total_reac_assoc = len(total_reac_assoc)
                try:
                    current_pwy_dict["completion_rate"][padmet_source].add(rxn_id)
                except KeyError:
                    current_pwy_dict["completion_rate"][padmet_source]= set([rxn_id])

            current_pwy_dict["total_reac_assoc"] = total_reac_assoc
            reac_assoc = [rlt.id_in for rlt in padmet.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway"]
            nb_reac_assoc = len(reac_assoc)
            current_rxn_dict["pathway_assoc"][pwy_id][padmet_source] = {"nb_total_rxn": nb_total_reac_assoc, "nb_found_rxn": nb_reac_assoc}
            current_pwy_dict["reac_assoc"][rxn_id].add(padmet_source)
                
        #get all genes associated to the reaction
        linked_rlt = [rlt for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]
        if linked_rlt:
            for rlt in linked_rlt:
                gene_id = rlt.id_out
                gene_node = padmet.dicOfNode[gene_id]
                current_padmet_dict["gene"][gene_id] = dict()

                if gene_id not in current_rxn_dict["gene_assoc"].keys():
                    current_rxn_dict["gene_assoc"][gene_id] = {padmet_source: dict()}
                if gene_id not in total_padmet_data["gene"].keys():
                    total_padmet_data["gene"][gene_id] = {"basic_attrib": dict(), "padmet_source":set(), "reac_assoc": dict(), "pathway_assoc":dict()}
                current_gene_dict = total_padmet_data["gene"][gene_id]
                current_gene_dict = update_basic_attrib(gene_node, current_gene_dict, padmet_source)      
                current_gene_dict["padmet_source"].add(padmet_source)
                
                
                for pwy_id, pwy_dict in current_rxn_dict["pathway_assoc"].items():
                    if padmet_source in pwy_dict.keys():
                        pwy_current_source = pwy_dict[padmet_source]
                        if pwy_id not in current_gene_dict["pathway_assoc"].keys():
                            current_gene_dict["pathway_assoc"][pwy_id] = {padmet_source: dict()}
                        current_gene_dict["pathway_assoc"][pwy_id][padmet_source] = dict(pwy_current_source)

                if rxn_id not in current_gene_dict["reac_assoc"].keys():
                    current_gene_dict["reac_assoc"][rxn_id] = {padmet_source: dict()}
                #a is_linked_to rlt have in misc a key "SOURCE:ASSIGNMENT"
                #the value can be only the source, ex: OUTPUT_PANTOGRAPH_X
                #or the source and the known assignment: SILI_ANNOTATION:EC-NUMBER
                sources = rlt.misc["SOURCE:ASSIGNMENT"]
                for source_data in sources:
                    #ValueError: can't split and get 2 value == no known assignment
                    try:
                        source, assignment = source_data.split(":")
                    except ValueError:
                        source = source_data
                        assignment = "n.a"
                    rec_data_node = next(node for node in list(padmet.dicOfNode.values()) if (node.type == "reconstructionData" and node.misc.get("SOURCE",["unknown-source"])[0] == source))
                    category = rec_data_node.misc.get("CATEGORY",["unknown-category"])[0].lower()
                    tool = rec_data_node.misc.get("TOOL",["unknown-tool"])[0].lower()
                    comment = rec_data_node.misc.get("COMMENT",["n.a"])[0].lower()
                    source = source.lower()
                    assignment = assignment.lower()
                    if re.findall("^output_.*_from_", source):
                        source = re.sub("^output_.*_from_", "", source)

                    #special patch for aureme workflow, if source start with 'output_*_"
                    dict_src_data = {"source": source, "assignment":assignment, "tool": tool, "comment": comment}
                    try:
                        current_rxn_dict["gene_assoc"][gene_id][padmet_source][category].append(dict_src_data)
                    except KeyError:
                        current_rxn_dict["gene_assoc"][gene_id][padmet_source][category] = list([dict_src_data])
                    try:
                        current_rxn_dict["padmet_gene"][padmet_source].add(gene_id)
                    except KeyError:
                        current_rxn_dict["padmet_gene"][padmet_source] = set([gene_id])

                    try:
                        current_gene_dict["reac_assoc"][rxn_id][padmet_source][category].append(dict_src_data)
                    except KeyError:
                        current_gene_dict["reac_assoc"][rxn_id][padmet_source][category] = list([dict_src_data])
                        
        reconstruction_data = [padmet.dicOfNode[rlt.id_out] for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"]
        #src_data = {category:{source:{comment:comment, tool:tool}}}
        
        for rec_data_node in reconstruction_data:
            #if found, lower to standardize
            category = rec_data_node.misc.get("CATEGORY",["unknown-category"])[0].lower()
            tool = rec_data_node.misc.get("TOOL",["unknown-tool"])[0].lower()
            if category == "manual" and tool == "unknown-tool":
                tool = "curation"
            comment = rec_data_node.misc.get("COMMENT",["n.a"])[0].lower()
            source = rec_data_node.misc.get("SOURCE",["unknown-source"])[0].lower()
            if re.findall("^output.*_from_", source):
                source = re.sub("^output_.*_from_", "", source)
            if category not in total_padmet_data["reconstruction"]["category"].keys():
                total_padmet_data["reconstruction"]["category"][category] = dict()
            if tool not in total_padmet_data["reconstruction"]["tool"].keys():
                total_padmet_data["reconstruction"]["tool"][tool] = dict()
            if comment not in total_padmet_data["reconstruction"]["comment"].keys():
                total_padmet_data["reconstruction"]["comment"][comment] = dict()
            if source not in total_padmet_data["reconstruction"]["source"].keys():
                total_padmet_data["reconstruction"]["source"][source] = dict()
            dict_src_data = {"category": category, "source": source, "tool": tool, "comment": comment}
            current_rxn_dict["reconstruction_data"][padmet_source].append(dict_src_data)

            try:
                total_padmet_data["reconstruction"]["category"][category][rxn_id].add(padmet_source)
            except KeyError:
                total_padmet_data["reconstruction"]["category"][category][rxn_id] = set([padmet_source])
            try:
                total_padmet_data["reconstruction"]["tool"][tool][rxn_id].add(padmet_source)
            except KeyError:
                total_padmet_data["reconstruction"]["tool"][tool][rxn_id] = set([padmet_source])
            try:
                total_padmet_data["reconstruction"]["comment"][comment][rxn_id].add(padmet_source)
            except KeyError:
                total_padmet_data["reconstruction"]["comment"][comment][rxn_id] = set([padmet_source])
            try:
                total_padmet_data["reconstruction"]["source"][source][rxn_id].add(padmet_source)
            except KeyError:
                total_padmet_data["reconstruction"]["source"][source][rxn_id] = set([padmet_source])

    return total_padmet_data

def update_basic_attrib(node, current_node_dict, padmet_source):
    """
    #TODO
    """
    for tag, list_of_v in node.misc.items():
        tag = tag.lower()
        if tag not in current_node_dict["basic_attrib"].keys():
            current_node_dict["basic_attrib"][tag] = dict()
        global_v = current_node_dict["basic_attrib"][tag]
        for v in list_of_v:
            v = v.lower()
            if tag == "molecular-weight": v = v.replace(" ","")
            try:
                global_v[v].add(padmet_source)
            except KeyError:
                global_v[v] = set([padmet_source])
    return current_node_dict

def reduce_padmet_data(total_padmet_data, verbose=False):
    """
    #TODO
    """
    all_padmet_source = set(total_padmet_data["padmet_source"].keys())
    if len(all_padmet_source) == 1:
        global_val = "_1_IN_1_"
    else:
        global_val = "_1_IN_ALL_"
    for rxn_id, rxn_data in total_padmet_data["reaction"].items():
        basic_attrib_dict = rxn_data["basic_attrib"]
        for tag, dict_val_padmet_source in basic_attrib_dict.items():
            for val, set_padmet_source in dict_val_padmet_source.items():
                if set_padmet_source == all_padmet_source:
                    dict_val_padmet_source[val] = global_val
        formula_dict = rxn_data["formula"]
        if len(formula_dict.keys()) == 1:
            formula_dict[list(formula_dict.keys())[0]] = global_val
            
                        

    for cpd_id, cpd_data in total_padmet_data["metabolite"].items():
        basic_attrib_dict = cpd_data["basic_attrib"]
        for tag, dict_val_padmet_source in basic_attrib_dict.items():
            for val, set_padmet_source in dict_val_padmet_source.items():
                if set_padmet_source == all_padmet_source:
                    dict_val_padmet_source[val] = global_val
        for key in ["consumed_by", "produced_by", "unknown_dir_by"]:
            set_padmet_source = cpd_data[key]
            if set_padmet_source == all_padmet_source:
                cpd_data[key] = global_val

    for gene_id, gene_data in total_padmet_data["gene"].items():
        basic_attrib_dict = gene_data["basic_attrib"]
        for tag, dict_val_padmet_source in basic_attrib_dict.items():
            for val, set_padmet_source in dict_val_padmet_source.items():
                if set_padmet_source == all_padmet_source:
                    dict_val_padmet_source[val] = global_val
        
    for pwy_id, pwy_data in total_padmet_data["pathway"].items():
        basic_attrib_dict = pwy_data["basic_attrib"]
        for tag, dict_val_padmet_source in basic_attrib_dict.items():
            for val, set_padmet_source in dict_val_padmet_source.items():
                if set_padmet_source == all_padmet_source:
                    dict_val_padmet_source[val] = global_val
        reac_assoc_dict = pwy_data["reac_assoc"]
        for rxn_id, set_padmet_source in reac_assoc_dict.items():
            if set_padmet_source == all_padmet_source:
                reac_assoc_dict[rxn_id] = global_val
        completion_rate_dict = pwy_data["completion_rate"]
        for padmet_source in all_padmet_source:
            try:
                rxn_assoc = completion_rate_dict[padmet_source]
            except KeyError:
                rxn_assoc = []

            try:
                rate = round(float(len(rxn_assoc))/float(len(pwy_data["total_reac_assoc"])),2)
            except ZeroDivisionError:
                rate = None
            completion_rate_dict[padmet_source] = rate
            
            

    return total_padmet_data
                    
    

def createDirectory(output, verbose=False):
    """
    create the folders genes, reactions, metabolites, pathways in the folder dirPath/
    if already exist, it will replace old folders (and delete old files)

    Parameters
    ----------
    output: str
        path to output folder
    """
    #simple check that dirPath is a dir:
    dirNames = ["genes","reactions","metabolites","pathways","organisms","navigation","files","misc"]
    #creatings the directory which will contains the wiki pages
    for d in dirNames:
        d_path = os.path.join(output,d)
        if not os.path.exists(d_path):
            if verbose: print("Creating directory: %s" %d_path)
            os.makedirs(d_path)
        else:
            if verbose: print("The directory %s already exist. Old pages will be deleted" %d_path)
            shutil.rmtree(d_path)
            os.makedirs(d_path)
    
    

def create_biological_page(category, page_id, page_dict_data, total_padmet_data, ext_link, output_file, padmetRef=None, verbose=False):
    """
    #TODO
    """
    total_padmet_source = list(total_padmet_data["padmet_source"].keys())
    to_collapse_cutoff = 20
    if verbose: print("%s: %s" %(category, page_id), end='\r')
    #stock in properties: all properties associated to the current page for semantic mediawiki
    #properties = [{{#set PROPERTY_X:VALUE_1|...|VALUE_N}}, ...]
    properties = []

    #ext_link is used to create external link to the database of reference if known
    if ext_link.get(category):
        dataInArray = ['[[Category:'+category+']]',
        '== '+category.capitalize()+' ['+ext_link.get(category)+page_id+' '+page_id+'] ==']
    else:
        dataInArray = ['[[Category:'+category+']]',
        '== '+category.capitalize()+' '+page_id+' ==']

    #basic_attrib contains data from misc dict of nodes ex for reaction: common-name, ec-number ...
    #page_dict_data["basic_attrib"] = {tag*: dict_val_padmet_source}
    #dict_va_padmet_source = {value: set_padmet_source}
    #ex: {'ec-number':{'ec-1.2.2: set([padmet_a, padmet_b])}}
    if "basic_attrib" in page_dict_data.keys():
        for tag, dict_val_padmet_source in page_dict_data["basic_attrib"].items():
            line = '* '+tag+':'
            dataInArray.append(line)
            #dont add smile as propertie because of invalid characters for SemanticMediaWiki
            if tag not in ["smiles"]:
                add_property(properties, tag, dict_val_padmet_source.keys())
            for val, set_padmet_source in dict_val_padmet_source.items():
                #for ec-number, remove 'ec-' and add external link to expasy
                if tag == "ec-number":
                    line = "** [http://enzyme.expasy.org/EC/"+val.replace("ec-","")+" "+val+"]"
                #for taxonomic-range from metacyc, add external link to biocyc
                elif tag == "taxonomic-range" and ext_link.get(category):
                    line = "** ["+ext_link.get(category)+val+" "+val+"]"
                #for inchi-key, remove "inchikey="
                elif tag == "inchi-key":
                    line = "** %s" %(val.replace("inchikey=",""))
                else:
                    line = "** %s" %val
                dataInArray.append(line)
                #when set_padmet_source == '_1_IN_ALL_': code means that the value is in all padmet/organisms
                if set_padmet_source == "_1_IN_ALL_":
                    line = "*** Value found in all organisms"
                    dataInArray.append(line)
                #when set_padmet_source != '_1_IN_ALL_' or '_1_IN_1': its a real set of padmet sources
                if set_padmet_source not in ["_1_IN_ALL_", "_1_IN_1_"]:
                    
                    line = "*** Value found in %s organism(s): %s" %(len(set_padmet_source), ";".join(["[[%s]]"%src for src in set_padmet_source]))
                    dataInArray.append(line)

    #For each category, extract in a specific way the information

    if category == "reaction":
        dataInArray.append("== Reaction formula ==")
        #add_property(properties, "formula", page_dict_data["formula"].keys())
        for formula, set_padmet_source in page_dict_data["formula"].items():
            line = "* "+formula
            dataInArray.append(line)
            if set_padmet_source == "_1_IN_ALL_":
                line = "** Same formula in all organisms"
                dataInArray.append(line)
            if set_padmet_source not in ["_1_IN_ALL_", "_1_IN_1_"]:
                line = "** Formula found in %s organism(s): %s" %(len(set_padmet_source), ";".join(set_padmet_source))
                dataInArray.append(line)

        if len(total_padmet_source) > 1:
            add_property(properties, "nb organism associated", [len(page_dict_data["padmet_source"])])
            min_gene_assoc = min([len(genes_set) for genes_set in page_dict_data["padmet_gene"].values()])
            max_gene_assoc = max([len(genes_set) for genes_set in page_dict_data["padmet_gene"].values()])
            avg_gene_assoc = round(sum([len(genes_set) for genes_set in page_dict_data["padmet_gene"].values()])/len(page_dict_data["padmet_gene"].values()),2)
            add_property(properties, "min gene associated", [min_gene_assoc])
            add_property(properties, "max gene associated", [max_gene_assoc])
            add_property(properties, "avg gene associated", [avg_gene_assoc])
            dataInArray.append('== Organism(s) associated with this reaction  ==')
            [dataInArray.append("* [["+padmet_source+"]]") for padmet_source in page_dict_data["padmet_source"]]
        else:
            add_property(properties, "nb gene associated", [len(page_dict_data["gene_assoc"].keys())])
            
        dataInArray.append('== Gene(s) associated with this reaction  ==')
        for gene_id, gene_dict in page_dict_data["gene_assoc"].items():
            line = "* Gene: [[%s]]" %gene_id
            dataInArray.append(line)
            if len(total_padmet_source) > 1:
                for padmet_source, category_dict in gene_dict.items():
                    line = "** In organism: [[%s]]" %padmet_source
                    dataInArray.append(line)
                    for category, category_data in category_dict.items():
                        line = "*** Category: [[%s]]" %category
                        dataInArray.append(line)
                        for data in category_data:
                            line = "**** Source: [[%s]], Tool: [[%s]], Assignment: %s, Comment: %s" %(data["source"], data["tool"], data["assignment"], data["comment"])
                            dataInArray.append(line)
            else:
                for category_dict in gene_dict.values():
                    for category, category_data in category_dict.items():
                        line = "** Category: [[%s]]" %category
                        dataInArray.append(line)
                        for data in category_data:
                            line = "*** Source: [[%s]], Tool: [[%s]], Assignment: %s, Comment: %s" %(data["source"], data["tool"], data["assignment"], data["comment"])
                            dataInArray.append(line)
        
        dataInArray.append('== Pathway(s)  ==')
        #set of pathways id associated to the reaction
        add_property(properties, "nb pathway associated", [len(page_dict_data["pathway_assoc"].keys())])
        for pwy_id, pwy_dict in page_dict_data["pathway_assoc"].items():
            if padmetRef and pwy_id in padmetRef.dicOfNode.keys():
                pwy_cname = padmetRef.dicOfNode[pwy_id].misc.get("COMMON-NAME",[None])[0]
                if pwy_cname:
                    line = "* [["+pwy_id+"]], "+pwy_cname+":"
            else:
                line = "* [["+pwy_id+"]]:"

            #if external link known, adding external link to the pathway
            #line: PWY_ID, common name, extern link to PWY
            if ext_link.get("pathway"):
                line += ' ['+ext_link.get("pathway")+pwy_id+' '+pwy_id+']'
            dataInArray.append(line)

            if len(total_padmet_source) > 1:
                for padmet_source, pwy_rxn in pwy_dict.items():
                    line = "** In organism: [[%s]]" %padmet_source
                    dataInArray.append(line)
                    line = "*** '''"+str(pwy_rxn["nb_found_rxn"])+"''' reactions found over '''"+str(pwy_rxn["nb_total_rxn"])+"''' reactions in the full pathway"
                    dataInArray.append(line)
            else:
                pwy_rxn = list(pwy_dict.values())[0]
                line = "** '''"+str(pwy_rxn["nb_found_rxn"])+"''' reactions found over '''"+str(pwy_rxn["nb_total_rxn"])+"''' reactions in the full pathway"
                dataInArray.append(line)

        dataInArray.append('== Reconstruction information  ==')
        if len(total_padmet_source) > 1:
            for tool, tool_dict in total_padmet_data["reconstruction"]["tool"].items():
                try:
                    tool_org_list = tool_dict[page_id]
                except KeyError:
                    tool_org_list = list()
                add_property(properties, "nb organism associated based on %s tool"%tool, [len(tool_org_list)])
            for category, category_dict in total_padmet_data["reconstruction"]["category"].items():
                try:
                    category_org_list = category_dict[page_id]
                except KeyError:
                    category_org_list = list()
                add_property(properties, "nb organism associated based on %s category"%category, [len(category_org_list)])
          
            for padmet_source, rec_data in page_dict_data["reconstruction_data"].items():
                line="* In organism: [[%s]]" %padmet_source
                dataInArray.append(line)

                for rec_dict in rec_data:
                    line = "** category: [[%s]]; source: [[%s]]; tool: [[%s]]; comment: %s" %(rec_dict["category"], rec_dict["source"], rec_dict["tool"], rec_dict["comment"])
                    dataInArray.append(line)
        else:
            rec_data = list(page_dict_data["reconstruction_data"].values())[0]
            add_property(properties, "reconstruction category", [rec_dict["category"] for rec_dict in rec_data])
            add_property(properties, "reconstruction tool", [rec_dict["tool"] for rec_dict in rec_data])
            add_property(properties, "reconstruction comment", [rec_dict["comment"] for rec_dict in rec_data])
            add_property(properties, "reconstruction source", [rec_dict["source"] for rec_dict in rec_data])
            for rec_dict in rec_data:
                line = "* category: [[%s]]; source: [[%s]]; tool: [[%s]]; comment: %s" %(rec_dict["category"], rec_dict["source"], rec_dict["tool"], rec_dict["comment"])
                dataInArray.append(line)

    elif category == "gene":
        dataInArray.append('== Organism(s) associated with this gene  ==')
        add_property(properties, "organism associated", page_dict_data["padmet_source"])
        for padmet_source in page_dict_data["padmet_source"]:
            line = "* [[%s]]" %padmet_source
            dataInArray.append(line)
        dataInArray.append("== Reaction(s) associated ==")
        add_property(properties, "nb reaction associated", [len(page_dict_data["reac_assoc"].keys())])
        for rxn_id, rxn_dict in page_dict_data["reac_assoc"].items():
            line = "* [[%s]]" %rxn_id
            dataInArray.append(line)
            if len(total_padmet_source) > 1:
                for padmet_source, link_data in rxn_dict.items():
                    line = "** In organism: [[%s]]" %padmet_source
                    dataInArray.append(line)
                    for category, category_dict in link_data.items():
                        line = "*** Category: [[%s]]" %category
                        dataInArray.append(line)
                        for source_dict in category_dict:
                            line = "**** source: [[%s]]; tool: [[%s]]; comment: %s" %(source_dict["source"], source_dict["tool"], source_dict["comment"])
                            dataInArray.append(line)
            else:
                link_data = list(rxn_dict.values())[0]
                for category, category_dict in link_data.items():
                    line = "** Category: [[%s]]" %category
                    dataInArray.append(line)
                    for source_dict in category_dict:
                        line = "*** source: [[%s]]; tool: [[%s]]; comment: %s" %(source_dict["source"], source_dict["tool"], source_dict["comment"])
                        dataInArray.append(line)

        if page_dict_data["pathway_assoc"]:
            add_property(properties, "nb pathway associated", [len(page_dict_data["pathway_assoc"].keys())])
            dataInArray.append("== Pathway(s) associated ==")
            for pwy_id, pwy_dict in page_dict_data["pathway_assoc"].items():
                line = "* [[%s]]" %pwy_id
                dataInArray.append(line)
                if len(total_padmet_source) > 1:
                    for padmet_source, pwy_rxn in pwy_dict.items():
                        line = "** In organism [[%s]]: '''%s''' reactions found over '''%s''' reactions in the full pathway" %(padmet_source, pwy_rxn["nb_found_rxn"], pwy_rxn["nb_total_rxn"])
                else:
                    pwy_rxn = list(pwy_dict.values())[0]
                    line = "** '''%s''' reactions found over '''%s''' reactions in the full pathway" %(pwy_rxn["nb_found_rxn"], pwy_rxn["nb_total_rxn"])
                dataInArray.append(line)

    elif category == "pathway":
        if len(total_padmet_source) > 1:
            dataInArray.append('== Organism(s) associated with this pathway  ==')
            add_property(properties,"nb organism associated", [len(page_dict_data["padmet_source"])])
            rate_list = [rate for rate in page_dict_data["completion_rate"].values()]
            try:
                min_completion_rate = min(rate_list)
                max_completion_rate = max(rate_list)
                avg_completion_rate = round(sum(rate_list)/len(rate_list),2)
            except TypeError:
                min_completion_rate = "n.a"
                max_completion_rate = "n.a"
                avg_completion_rate = "n.a"
            add_property(properties,"min completion rate", [min_completion_rate])
            add_property(properties,"max completion rate", [max_completion_rate])
            add_property(properties,"avg completion rate", [avg_completion_rate])
            for padmet_source in page_dict_data["padmet_source"]:
                line = "* [[%s]]" %padmet_source
                dataInArray.append(line)
        else:
            nb_reaction_found = len(page_dict_data["reac_assoc"].keys())
            completion_rate = page_dict_data["completion_rate"][total_padmet_source[0]]
            if completion_rate is None:
                completion_rate = "n.a"
            add_property(properties,"nb reaction found", [nb_reaction_found])
            add_property(properties,"completion rate", [completion_rate])


        dataInArray.append("== Reaction(s) found ==")
        for rxn_id, set_padmet_source in page_dict_data["reac_assoc"].items():
            line = "* [[%s]]" %rxn_id
            dataInArray.append(line)
            if set_padmet_source == "_1_IN_ALL_":
                line = "** Reaction found in all organisms"
                dataInArray.append(line)
            if set_padmet_source not in ["_1_IN_ALL_", "_1_IN_1_"]:
                line = "** Reaction found in %s organism(s): %s" %(len(set_padmet_source), ";".join(set_padmet_source))
                dataInArray.append(line)
        total_rxn_assoc = page_dict_data["total_reac_assoc"]
        if len(total_padmet_source) > 1:
            dataInArray.append("== Reaction(s) not found in any organism ==")
        else:
            dataInArray.append("== Reaction(s) not found ==")
        if total_rxn_assoc:
            add_property(properties,"nb total reaction", [len(total_rxn_assoc)])
            rxn_not_found = page_dict_data["total_reac_assoc"].difference(page_dict_data["reac_assoc"].keys())
            if rxn_not_found:
                for rxn_id in rxn_not_found:
                    line = "* [{0}{1} {1}]".format(ext_link.get(category), rxn_id)
                    dataInArray.append(line)
            else:
                line = "All reactions of this pathways are in present"
                dataInArray.append(line)
        else:
            add_property(properties,"nb total reaction", ["n.a"])
            line = "No padmetRef was given during wikipage creation or pathway not in metacyc, data not available"
            dataInArray.append(line)

    elif category == "metabolite":
        dataInArray.append("== Reaction(s) known to consume the compound ==")
        if page_dict_data["consumed_by"]:
            #add_property(properties, "consumed by", page_dict_data["consumed_by"].keys())
            for rxn_id, set_padmet_source in page_dict_data["consumed_by"].items():
                line = "* [[%s]]" %rxn_id
                dataInArray.append(line)
                if len(total_padmet_source) > 1:
                    for padmet_source in set_padmet_source:
                        line = "** In organism [[%s]]" %padmet_source
                        dataInArray.append(line)
                        
        dataInArray.append("== Reaction(s) known to produce the compound ==")
        if page_dict_data["produced_by"]:
            #add_property(properties, "produced by", page_dict_data["produced_by"].keys())
            for rxn_id, set_padmet_source in page_dict_data["produced_by"].items():
                line = "* [[%s]]" %rxn_id
                dataInArray.append(line)
                if len(total_padmet_source) > 1:
                    for padmet_source in set_padmet_source:
                        line = "** In organism [[%s]]" %padmet_source
                        dataInArray.append(line)

        dataInArray.append("== Reaction(s) of unknown directionality ==")
        if page_dict_data["unknown_dir_by"]:
            #add_property(properties, "associated to unknown reaction direction", page_dict_data["unknown_dir_by"].keys())
            for rxn_id, set_padmet_source in page_dict_data["unknown_dir_by"].items():
                line = "* [[%s]]" %rxn_id
                dataInArray.append(line)
                if len(total_padmet_source) > 1:
                    for padmet_source in set_padmet_source:
                        line = "** In organism [[%s]]" %padmet_source
                        dataInArray.append(line)

    elif category == "organism":
        dataInArray.append("== Summary ==")
        nb_rxn, nb_gene, nb_pwy = [len(page_dict_data[i].keys()) for i in ["reaction","gene","pathway"]]
        add_property(properties, "nb reaction associated", [nb_rxn])
        add_property(properties, "nb gene associated", [nb_gene])
        add_property(properties, "nb pathway associated", [nb_pwy])
        dataInArray.append("* %s reaction(s)" %nb_rxn)
        dataInArray.append("* %s gene(s)" %nb_gene)
        dataInArray.append("* %s pathway(s)" %nb_pwy)
        dataInArray.append("== Reaction(s) associated ==")
        for rxn_id in page_dict_data["reaction"].keys():
            line = "* [[%s]]" %rxn_id
            dataInArray.append(line)
        dataInArray.append("== Gene(s) associated ==")
        for gene_id in page_dict_data["gene"].keys():
            line = "* [[%s]]" %gene_id
            dataInArray.append(line)
        dataInArray.append("== Pathway(s) associated ==")
        for pwy_id in page_dict_data["pathway"].keys():
            line = "* [[%s]]" %pwy_id
            dataInArray.append(line)
        
    """
    elif category == "diffCondition":
        genes_up = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_de_in" and float(rlt.misc["LOG2FC"][0]) > float(0)])
        genes_down = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_de_in" and float(rlt.misc["LOG2FC"][0]) < float(0)])
        rxns_up = set()
        rxns_down = set()
        rxns_unknown = set()
        pwys_up = set()
        pwys_down = set()
        pwys_unknown = set()
        dataInArray.append("== %s Gene(s) up-regulated ==" %len(genes_up))
        for gene_id in genes_up:
            dataInArray.append("* [[%s]]" %gene_id)
            rxns = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut[gene_id] if rlt.type == "is_linked_to"])
            rxns_up.update(rxns)
        dataInArray.append("== %s Gene(s) down-regulated ==" %len(genes_down))
        for gene_id in genes_down:
            dataInArray.append("* [[%s]]" %gene_id)
            rxns = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut[gene_id] if rlt.type == "is_linked_to"])
            rxns_down.update(rxns)

        rxns_unknown = rxns_up.intersection(rxns_down)
        if rxns_unknown:
            rxns_up.difference_update(rxns_down)
            rxns_down.difference_update(rxns_up)

        dataInArray.append("== %s Reaction(s) exclusively linked to up-regulated gene(s) ==" %len(rxns_up))
        for rxn_id in rxns_up:
            dataInArray.append("* [[%s]]" %rxn_id)
            pwys = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway" and rlt.id_out in total_pwy_id])
            pwys_up.update(pwys)
        dataInArray.append("== %s Reaction(s) exclusively linked to down-regulated gene(s) ==" %len(rxns_down))
        for rxn_id in rxns_down:
            dataInArray.append("* [[%s]]" %rxn_id)
            pwys = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway" and rlt.id_out in total_pwy_id])
            pwys_down.update(pwys)
        dataInArray.append("== %s Reaction(s) linked to up and down regulated genes ==" %len(rxns_unknown))
        for rxn_id in rxns_unknown:
            dataInArray.append("* [[%s]]" %rxn_id)
            pwys = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway" and rlt.id_out in total_pwy_id])
            pwys_unknown.update(pwys)

        pwys_unknown.update(pwys_up.intersection(pwys_down))
        if pwys_unknown:
            pwys_up.difference_update(pwys_down)
            pwys_down.difference_update(pwys_up)

        dataInArray.append("== %s Pathway(s) only with reactions exclusively linked to up-regulated gene(s) ==" %len(pwys_up))
        for pwy_id in pwys_up:
            dataInArray.append("* [[%s]]" %pwy_id)
        dataInArray.append("== %s Pathway(s) only with reactions exclusively linked to down-regulated gene(s) ==" %len(pwys_down))
        for pwy_id in pwys_down:
            dataInArray.append("* [[%s]]" %pwy_id)
        dataInArray.append("== %s Pathway(s) with reactions linked to up and donw regulated genes ==" %len(pwys_unknown))
        for pwy_id in pwys_unknown:
            dataInArray.append("* [[%s]]" %pwy_id)

        add_property(properties, "gene up-regulated", genes_up)
        add_property(properties, "gene down-regulated", genes_down)
        add_property(properties, "reaction up-regulated", rxns_up)
        add_property(properties, "reaction down-regulated", rxns_down)
        add_property(properties, "reaction up-and-down-regulated", rxns_down)
        add_property(properties, "pathway up-regulated", pwys_up)
        add_property(properties, "pathway down-regulated", pwys_down)
        add_property(properties, "pathway up-and-down-regulated", pwys_unknown)
    """
    if "xref" in page_dict_data.keys():
        dataInArray.append('== External links  ==')
        for db, ids in page_dict_data["xref"].items():
            xrefLink(dataInArray, db, ids)

    title_indices = [i for i, x in enumerate(dataInArray) if str(x).startswith("==")]
    for i in range(len(title_indices)):
        try:
            if title_indices[i+1] - title_indices[i] > to_collapse_cutoff:
                dataInArray[title_indices[i]+1: title_indices[i+1]] = add_collapsible(dataInArray[title_indices[i]+1: title_indices[i+1]])
                title_indices = [i for i, x in enumerate(dataInArray) if str(x).startswith("==")]
        except IndexError:
            if len(dataInArray) - title_indices[i] > to_collapse_cutoff:
                dataInArray[title_indices[i]+1:] = add_collapsible(dataInArray[title_indices[i]+1:])
        

    dataInArray.extend(properties)
    with open(output_file,'w', encoding="utf8") as f_handler:
        for line in dataInArray:
            f_handler.write(line+"\n")
    """                
    for i in dataInArray:
        print(i)
    print("\n")
    """
def add_property(properties, prop_id, prop_values):
    """
    #TODO
    """
    prop_values = set([str(i) for i in prop_values])
    start_line = "{{#set: "+prop_id+"="
    values_part = "|".join(prop_values)
    end_line = "}}"
    toInsert = start_line + values_part + end_line
    properties.append(toInsert)

def add_collapsible(text_array, title=None):
    """
    #TODO
    """
    collapsible = ['<div class="toccolours mw-collapsible mw-collapsed" style="width:100%; overflow:auto;">']
    if title:
        collapsible.append('<div style="font-weight:bold;line-height:1.6;">%s</div>' %title)
        collapsible.append('<div class="mw-collapsible-content">')
        collapsible.extend(text_array)
        collapsible.append('</div></div>')
    else:
        collapsible.extend(text_array)
        collapsible.append('</div>')
    return collapsible        

def xrefLink(dataInArray, db, ids):
    """
    #TODO
    """
    if db == "METACYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)
    elif db == "UNIPROT":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.uniprot.org/uniprot/"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "KEGG" or db.startswith("LIGAND"):
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.genome.jp/dbget-bin/www_bget?"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "RHEA":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.ebi.ac.uk/rhea/reaction.xhtml?id="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "WIKIPEDIA":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://en.wikipedia.org/wiki/"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "CHEBI":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.ebi.ac.uk/chebi/searchId.do?chebiId="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "PUBCHEM":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "ECOCYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://metacyc.org/ECOLI/NEW-IMAGE?object="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "CHEMSPIDER":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.chemspider.com/Chemical-Structure."+_id+".html "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "umbbd-compounds":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://umbbd.ethz.ch/servlets/pageservlet?ptype=c&compID="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "ARACYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://metacyc.org/ARA/NEW-IMAGE?object="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "PIR":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pir.georgetown.edu/cgi-bin/nbrfget?uid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "NCI":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://cactus.nci.nih.gov/ncidb2.2/?nsc="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "knapsacK":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://kanaya.naist.jp/knapsack_jsp/information.jsp?word="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    else:
        for _id in ids:
            toInsert = "* "+db+" : "+_id
            dataInArray.append(toInsert)


def create_navigation_page(total_padmet_data, navigation_folder, verbose=False):
    """
    #TODO
    """
    total_padmet_source = total_padmet_data["padmet_source"].keys()
    sideBarData = ["* navigation",
                   "** mainpage|mainpage-description",
                   "** Special:Ask|SMW-Ask",
                   "** workflow|workflow command history",
                   "** randompage-url|randompage",
                   "** Special:ListFiles|Files",
                   "* Metabolic network components"]

    category = "organism"
    sideBarData.append("** Category:{0}|{0}".format(category))
    fileName = os.path.join(navigation_folder, "Category:%s" %category)
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:%s]]" %category,
                   "| ?nb reaction associated",
                   "| ?nb gene associated", 
                   "| ?nb pathway associated",
                   "|sort=nb reaction associated", 
                   "|order=descending",
                   "}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "reaction"
    sideBarData.append("** Category:{0}|{0}".format(category))
    fileName = os.path.join(navigation_folder, "Category:%s" %category)
    if verbose: print("Category: %s" %category)
    
    if len(total_padmet_source) > 1:
        dataInArray = ["{{#ask: [[Category:%s]]"%category,
                       "| ?common-name",
                       "| ?ec-number",
                       "| ?nb organism associated",
                       "| ?min gene associated",
                       "| ?max gene associated",
                       "| ?avg gene associated",
                       "|sort=nb organism associated", 
                       "|order=descending"]

        dataInArray.extend(["| ?nb organism associated based on %s tool"%tool for tool in total_padmet_data["reconstruction"]["tool"].keys()]+["}}"])
    else:
        dataInArray = ["{{#ask: [[Category:%s]]"%category,
                       "| ?common-name",
                       "| ?ec-number",
                       "| ?nb gene associated",
                       "| ?nb pathway associated",
                       "| ?nb reconstruction source",
                       "| ?reconstruction category",
                       "| ?reconstruction tool",
                       "| ?reconstruction comment",
                       "}}"]
        
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "gene"
    sideBarData.append("** Category:{0}|{0}".format(category))
    fileName = os.path.join(navigation_folder, "Category:%s" %category)
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:%s]]"%category,
                   "| ?organism associated",
                   "| ?nb reaction associated",
                   "| ?nb pathway associated",
                   "|sort=nb reaction associated", 
                   "|order=descending",
                   "}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "pathway"
    sideBarData.append("** Category:{0}|{0}".format(category))
    fileName = os.path.join(navigation_folder, "Category:%s" %category)
    if verbose: print("Category: %s" %category)
    if len(total_padmet_source) > 1:
        dataInArray = ["{{#ask: [[Category:%s]]"%category,
                       "| ?common-name",
                       "| ?nb organism associated",
                       "| ?nb total reaction",
                       "| ?min completion rate",
                       "| ?max completion rate",
                       "| ?avg completion rate",
                       "|sort=avg completion rate", 
                       "|order=descending",
                       "}}"]
    else:
        dataInArray = ["{{#ask: [[Category:%s]]"%category,
                       "| ?common-name",
                       "| ?nb reaction found",
                       "| ?nb total reaction",
                       "| ?completion rate",
                       "|sort=completion rate, nb total reaction", 
                       "|order=descending",
                       "}}"]
        
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "metabolite"
    sideBarData.append("** Category:{0}|{0}".format(category))
    fileName = os.path.join(navigation_folder, "Category:%s" %category)
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:%s]]"%category,
                   "| ?common-name",
                   "}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")
    """
    if all_diffCdt:
        category = "diffCondition"
        sideBarData.append("* Transcriptomic data")
        sideBarData.append("** Category:"+category+"|"+category)
        fileName = output_folder+"Category:diffCondition"
        if verbose: print("Category: %s" %category)
        dataInArray = ["{{#ask: [[Category:diffCondition]]","| ?gene up-regulated","| ?reaction up-regulated","| ?pathway up-regulated","| ?gene down-regulated","| ?reaction down-regulated","| ?pathway down-regulated","| ?reaction up-and-down-regulated","| ?pathway up-and-down-regulated""}}"]
        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")
    """
    sideBarData.append("* Reconstruction categories")
    all_categories = total_padmet_data["reconstruction"]["category"].keys()
    [sideBarData.append("** "+rec_category+"|"+rec_category) for rec_category in sorted(all_categories)]
    for rec_category in all_categories:
        fileName = os.path.join(navigation_folder, rec_category)
        if verbose: print("Reconstruction category: %s" %rec_category)
        if len(total_padmet_source) > 1:
            dataInArray = ["{{#ask: [[Category:reaction]]",
                           "[[nb organism associated based on %s category::>1]]"%rec_category,
                           "| ?common-name",
                           "| ?ec-number",
                           "| ?nb organism associated",
                           "| ?min gene associated",
                           "| ?max gene associated",
                           "| ?avg gene associated",
                           "| ?nb organism associated based on %s category"%rec_category, 
                           "|sort=nb organism associated based on %s category"%rec_category, 
                           "|order=descending",
                           "}}"]
        else:
            dataInArray = ["{{#ask: [[Category:reaction]] [[reconstruction category::%s]]" %rec_category,
                           "| ?common-name",
                           "| ?ec-number",
                           "| ?reconstruction tool",
                           "| ?reconstruction source",
                           "| ?reconstruction comment",
                           "| ?nb gene associated",
                           "| ?nb pathway associated",
                           "}}"]

        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")

    sideBarData.append("* Reconstruction tools")
    all_tools = total_padmet_data["reconstruction"]["tool"].keys()
    [sideBarData.append("** "+rec_tool+"|"+rec_tool) for rec_tool in sorted(all_tools)]
    for rec_tool in all_tools:
        fileName = os.path.join(navigation_folder, rec_tool)
        if verbose: print("Reconstruction tool: %s" %rec_tool)
        if len(total_padmet_source) > 1:
            dataInArray = ["{{#ask: [[Category:reaction]]",
                           "[[nb organism associated based on %s tool::>1]]"%rec_tool,
                           "| ?common-name",
                           "| ?ec-number",
                           "| ?nb organism associated",
                           "| ?min gene associated",
                           "| ?max gene associated",
                           "| ?avg gene associated",
                           "| ?nb organism associated based on %s tool"%rec_tool, 
                           "|sort=nb organism associated based on %s tool"%rec_tool, 
                           "|order=descending",
                           "}}"]
        else:
            dataInArray = ["{{#ask: [[Category:reaction]] [[reconstruction tool::%s]]" %rec_tool,
                           "| ?common-name",
                           "| ?ec-number",
                           "| ?reconstruction category",
                           "| ?reconstruction source",
                           "| ?reconstruction comment",
                           "| ?nb gene associated",
                           "| ?nb pathway associated",
                           "}}"]

        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")

    """
    sideBarData.append("* Reconstruction sources")
    all_sources = total_padmet_data["reconstruction"]["source"].keys()
    [sideBarData.append("** "+rec_source+"|"+rec_source) for rec_source in sorted(all_sources)]

    for rec_source in set(all_sources).difference(total_padmet_data["padmet_source"].keys()):
        fileName = os.path.join(navigation_folder, rec_source)
        if verbose: print("Reconstruction source: %s" %rec_source)
        dataInArray = ["{{#ask: [[Category:Reaction]] [[reconstruction source::"+rec_source+"]]","| ?COMMON NAME","| ?ec number",
                          "| ?reconstruction category","| ?reconstruction tool","| ?reconstruction source","| ?reconstruction comment","| ?gene associated","| ?in pathway","}}"]
        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")
    """
    if verbose: print("SideBar page")
    fileName = os.path.join(navigation_folder, "MediaWiki:Sidebar")
    with open(fileName, 'w') as f:
        for line in sideBarData:
            f.write(line+"\n")


def create_venn(total_padmet_data, output_file, verbose=False):
    """
    #TODO
    """
    if verbose: print("Venn Diagramm")
               
    categories = ["annotation", "orthology", "manual", "gap-filling"]
    reactions_list = [total_padmet_data["reconstruction"]["category"].get(cat, {}).keys() for cat in categories]
    labels = get_labels(reactions_list)
    

    #print("venn debug")
    #print(labels)
    #print(categories)

    fig, ax = venn4(labels, categories)
    fig.savefig(output_file)

def copy_io_files():
    """
    """
    #toDo in futur
        
def create_main(total_padmet_data, wiki_id, output_file, verbose=False):
    """
    #TODO
    """
    if verbose: print("Main page")
    ### create main page
    all_rxns = len(total_padmet_data["reaction"].keys())
    all_cpds = len(total_padmet_data["metabolite"].keys())
    all_pwys = len(total_padmet_data["pathway"].keys())
    all_genes = len(total_padmet_data["gene"].keys())

    for line in main_template:
        main_template[main_template.index(line)] = line.replace("MODEL_ID",wiki_id)
    final_network_index = main_template.index([line for line in main_template if line.startswith("The automatic")][0])
    main_template[final_network_index] = main_template[final_network_index].replace("NB_RXN", str(all_rxns)).replace("NB_CPD", str(all_cpds)).replace("NB_PWY", str(all_pwys)).replace("NB_GENE", str(all_genes))

    """
    index += 1

    index += 1        

    """
    index = 1
    for category, rxns in total_padmet_data["reconstruction"]["category"].items():
        if category == "annotation":
            main_template.insert(final_network_index+index, "* Based on annotation data:")
            index += 1
            main_template.insert(final_network_index+index, "** Creation of a metabolic network containing %s reactions" %len(rxns))
            index += 1
        elif category == "orthology":
            main_template.insert(final_network_index+index, "* Based on orthology data:")
            index += 1
            main_template.insert(final_network_index+index, "** Creation of a global metabolic network containing %s reactions" %len(rxns))
            index += 1
            for src,v in total_padmet_data["reconstruction"]["source"].items():
                if src.startswith("output_pantograph_") or src.startswith("output_orthofinder_"):
                    src = src.replace("output_pantograph_","").replace("output_orthofinder_","")
                    main_template.insert(final_network_index+index, "*** From template ''"+src+"'' creation of a metabolic network containing: %s reactions"%len(v.keys()))
                    index += 1
        elif category == "manual":
            main_template.insert(final_network_index+index, "* Based on expertise:")
            index += 1
            main_template.insert(final_network_index+index, "** %s reaction(s) added"%len(rxns))
            index += 1
        elif category == "gap-filling":
            main_template.insert(final_network_index+index, "* Based on gap-filling:")
            index += 1
            main_template.insert(final_network_index+index, "*** %s reaction(s) added"%len(rxns))
            index += 1
    main_template.insert(final_network_index+index, "")
    index += 1
    main_template.insert(final_network_index+index, "*List of tools used:")
    index += 1
    for tool in total_padmet_data["reconstruction"]["tool"].keys():
        if tool == "pantograph":
            main_template.insert(final_network_index+index, "** Tool: [http://pathtastic.gforge.inria.fr Pantograph]")
        elif tool == "pathwaytools":
            main_template.insert(final_network_index+index, "** Tool: [http://bioinformatics.ai.sri.com/ptools/ PathwayTools]")
        elif tool == "meneco":
            main_template.insert(final_network_index+index, "** Tool: [https://pypi.python.org/pypi/meneco Meneco]")
        elif tool == "orthofinder":
            main_template.insert(final_network_index+index, "** Tool: [https://github.com/davidemms/OrthoFinder Orthofinder]")
        index += 1
        
            

    with open(output_file,'w') as f:
        for line in main_template:
            f.write(line+"\n")                


def get_labels(data, fill=["number"]):
    """    
    get a dict of labels for groups in data
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
    Out[12]:
    {'001': '0',
    '010': '5',
    '011': '0',
    '100': '3',
    '101': '2',
    '110': '2',
    '111': '3'}
    
    Parameters
    ----------    
    data: list
        data to get label for
    fill: 
        ["number"|"logic"|"percent"]

    Returns
    -------
    dict:
        a dict of labels for different sets
    """
    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        data_size = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)

    return labels

def venn4(labels, names=['A', 'B', 'C', 'D'], **options):
    """
    plots a 4-set Venn diagram

    Parameters
    ----------    
    labels: dict
        a label dict where keys are identified via binary codes ('0001', '0010', '0100', ...),
        hence a valid set could look like: {'0001': 'text 1', '0010': 'text 2', '0100': 'text 3', ...}.
        unmentioned codes are considered as ''.
    names: list
        group names

    Returns
    -------
    set
        (Figure, AxesSubplot), pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(4)])
    figsize = options.get('figsize', (12, 12))
    dpi = options.get('dpi', 96)
    
    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)
    
    # body   
    draw_ellipse(fig, ax, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
    draw_ellipse(fig, ax, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
    draw_ellipse(fig, ax, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
    draw_ellipse(fig, ax, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
    draw_text(fig, ax, 0.85, 0.42, labels.get('0001', ''))
    draw_text(fig, ax, 0.68, 0.72, labels.get('0010', ''))
    draw_text(fig, ax, 0.77, 0.59, labels.get('0011', ''))
    draw_text(fig, ax, 0.32, 0.72, labels.get('0100', ''))
    draw_text(fig, ax, 0.71, 0.30, labels.get('0101', ''))
    draw_text(fig, ax, 0.50, 0.66, labels.get('0110', ''))
    draw_text(fig, ax, 0.65, 0.50, labels.get('0111', ''))
    draw_text(fig, ax, 0.14, 0.42, labels.get('1000', ''))
    draw_text(fig, ax, 0.50, 0.17, labels.get('1001', ''))
    draw_text(fig, ax, 0.29, 0.30, labels.get('1010', ''))
    draw_text(fig, ax, 0.39, 0.24, labels.get('1011', ''))
    draw_text(fig, ax, 0.23, 0.59, labels.get('1100', ''))
    draw_text(fig, ax, 0.61, 0.24, labels.get('1101', ''))
    draw_text(fig, ax, 0.35, 0.50, labels.get('1110', ''))
    draw_text(fig, ax, 0.50, 0.38, labels.get('1111', ''))
    
    # legend
    draw_text(fig, ax, 0.13, 0.18, names[0], colors[0])
    draw_text(fig, ax, 0.18, 0.83, names[1], colors[1])
    draw_text(fig, ax, 0.82, 0.83, names[2], colors[2])
    draw_text(fig, ax, 0.87, 0.18, names[3], colors[3])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    return fig, ax

def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)

def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1]):
    ax.text(
        x, y, text,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=14,
        color=color)

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


main_template = ["== MODEL_ID description ==",
     "== Automatic reconstruction with [http://aureme.genouest.org AuReMe] ==",
     "Model summary: [[MEDIA:summary.txt|summary]]",
     "",
     "Download '''AuReMe''' Input/Output [LINK OR MEDIA data]",
     "",
     "The automatic reconstruction of ''MODEL_NAME'' results to a Genome scale [[MEDIA:model.xml|Model]] containing NB_RXN reactions, NB_CPD metabolites, NB_GENE genes and NB_PWY pathways. This GeM was obtained based on the following sources:",
     "",
     "[[FILE:venn.png|frameless|border]]",
     "",
     "== Collaborative curation == ",
     "* Suggest reactions to add or remove:",
     "** Download this [[MEDIA:Add_delete_reaction.csv|form]]",
     "* Suggest new reactions to create and add:",
     "** Download this [[MEDIA:Reaction_creator.csv|form]]",
     "* '''Follow the examples given in the form(s) to correctly share your suggestions'''",
     "* Send the filled form(s) to: CONTACT_MAIL"]

def create_log_page(log_file, output_folder):
    """
    #TODO
    """
    cmd_regex = '--cmd=\"(.*)\"'
    fileName = os.path.join(output_folder,"workflow")
    log_page = ["=Workflow command history=","","==Command sequence=="]
    with open(log_file, 'r') as f:
        log_data = [line for line in f.read().splitlines() if not line.startswith("#")]
    for cmd_line in log_data:
        re_result = re.search(cmd_regex, cmd_line)
        if re_result:
            full_cmd = re_result.groups(1)[0]
            cmd = full_cmd.split(" ")[0]
            cmd_label,desc = get_cmd_label(cmd)
            if cmd_label and desc:
                if cmd == "curation":
                    curation_file = re.search('DATA=(.*\.\w+)\s?',full_cmd).groups(1)[0]
                    desc = desc.replace('FORM_FILE_NAME', curation_file)
                log_page.append("* '''%s''':" %cmd_label)
                log_page.append("''%s''" %desc)
    log_page.extend(["==Downloads==","You can download the [[MEDIA:log.txt|command log file here]]"])
    #todo downalosalsjka
    with open(fileName, 'w') as f:
        for line in log_page:
            f.write(line+"\n")
        
                 

def get_cmd_label(cmd):
    """
    #TODO
    """
    cmd_label_dict = {'getdb':{'CMD_LABEL':'Get database','DESC':"Display the available reference databases"},
                      'check_input':{'CMD_LABEL':'Check input','DESC':"Check the validity, consistency and presence of input files"},
                      'check_studied_organism_input':{'CMD_LABEL':'Check studied organism input','DESC':"Check if FAA or GBK was given for the studied organism"},
                      'check_model_organism_input':{'CMD_LABEL':'Check model organism input','DESC':"Check (if existing) each folder in orthology based reconstruction"},
                      'sbml_validity':{'CMD_LABEL':"Check SBML validity", 'DESC':"Check the SBML validity of $(METABOLIC_MODEL)."},
                      'faa_validity':{'CMD_LABEL':"Check Fasta validity", 'DESC':"Check if genes IDs in the model metabolic network are the same than in the Fasta file.<BR>If the rate of validated genes IDs is lower than the CUTOFF (see config.txt), raise an error and break the workflow."},
                      'check_gap_filling_input':{'CMD_LABEL':"Check inputs for gap-filling", 'DESC':"Check the seeds, targets and artefacts files."},
                      'curation':{'CMD_LABEL':"Manual curation",'DESC':"Apply the curation described in the form file FORM_FILE_NAME."},
                      'annotation_based':{'CMD_LABEL':"Annotation based reconstruction", 'DESC':"Extract network data from Pathway Tools annotation output."},
                      'pathwaytools':{'CMD_LABEL':"Pathway Tools", 'DESC':"For each folder in /annotation_based_reconstruction, check or generate from pgdb files (.dat) the SBML file."},
                      'orthology_based':{'CMD_LABEL':"Orthology based reconstruction", 'DESC':"Run the orthology based reconstruction."},
                      'inparanoid':{'CMD_LABEL':"Inparanoid", 'DESC':"Run Inparanoid : search for orthologs."},
                      'OMCL':{'CMD_LABEL':"OrthoMCL", 'DESC':"Run OrthoMCL : search for orthologs."},
                      'mp_pantograph' :{'CMD_LABEL':"Pantograph multiprocessed", 'DESC':"Run Pantograph : merges OrthoMCL and inparanoid results.<BR>Quickened by mulitprocessing."},
                      'pantograph':{'CMD_LABEL':"Pantograph", 'DESC':"Run Pantograph : merges OrthoMCL and inparanoid results."},
                      'draft':{'CMD_LABEL':"Create draft network", 'DESC':"Merges all available networks from the /networks directory into one metabolic network.<BR>Merge all data on the studied species."},
                      'gap_filling':{'CMD_LABEL':"Run gap-filling", 'DESC':"Calculate the gap-filling solution and generate the metabolic network, completed with the gap-filling solution."},
                      'gap_filling_solution':{'CMD_LABEL':"Calculate gap-filling solution", 'DESC':"Only calculate the gap-filling solution."},
                      'meneco':{'CMD_LABEL':"Meneco", 'DESC':"Run Meneco : a gap-filling reconstruction method."},
                      'final':{'CMD_LABEL':"Final network", 'DESC':"Generate the final metabolic network, once applyed all the reconstruction methods."},
                      'gbk_to_faa':{'CMD_LABEL':"GBK to Fasta", 'DESC':"Export a GeneBank (.gbk) file in Fasta (faa) format."},
                      'pgdb_to_padmet':{'CMD_LABEL':"PGDB to PADMet", 'DESC':"Export a PGDB (.dat Pathway Tools files) in PADMet (.padmet) format."},
                      'padmet_to_sbml':{'CMD_LABEL':"PADMet to SBML", 'DESC':"Export a PADMet (.padmet) file in the SBML format."},
                      'compounds_to_sbml':{'CMD_LABEL':"Compounds to SBML", 'DESC':"Export a list of compounds (.txt) in the SBML format."},
                      'sbml_mapping':{'CMD_LABEL':"SBML mapping", 'DESC':"Map an SBML file (all the entities IDs) to a reference database (specified in config.txt)."},
                      'get_medium':{'CMD_LABEL':"Get medium", 'DESC':"Display the studied species specified growth medium."},
                      'set_medium':{'CMD_LABEL':"Set medium", 'DESC':"Set the growth medium for the studied species."},
                      'del_medium':{'CMD_LABEL':"Delete medium", 'DESC':"Delete the growth medium for the studied species."},
                      'get_compart':{'CMD_LABEL':"Get compartments", 'DESC':"Display all the compartments of the metabolic network."},
                      'del_compart':{'CMD_LABEL':"Delete compartment", 'DESC':"Remove a compartment from the metabolic network."},
                      'change_compart':{'CMD_LABEL':"Change compartment", 'DESC':"Modify a compartment in the metabolic network."},
                      'report':{'CMD_LABEL':"Report", 'DESC':"Generate reports on the metabolic network reconstruction."},
                      'set_fba':{'CMD_LABEL':"Set FBA", 'DESC':"Set the biomass reaction to run flux balance analysis on the network."},
                      'summary':{'CMD_LABEL':"Test FBA", 'DESC':"Run flux balance analysis on the network."},
                      'menecheck':{'CMD_LABEL':"Menecheck", 'DESC':"Run topological analysis on the network."},
                      'shogen':{'CMD_LABEL':"Shogen", 'DESC':"Run Shogen : find shortest genome segments that regulate metabolic pathways."},
                      'wiki_pages':{'CMD_LABEL':"Create Wiki pages", 'DESC':"Create Wiki pages to display the metabolic network."},
                      'wiki_run':{'CMD_LABEL':"Run Wiki", 'DESC':"Create a Docker container for the Wiki."},
                      'wiki_init':{'CMD_LABEL':"Wiki initialization", 'DESC':"Send data on the metabolic network to the Docker container to fill in the Wiki."},
                      'send_all_page':{'CMD_LABEL':"Send all pages to Wiki", 'DESC':"Send all the generated pages on the metabolic network to the Wiki."},
                      'tsv':{'CMD_LABEL':"PADMet to tsv",'DESC':"Convert a PADMet (.padmet) file to a tsv files (for Askomics)."}
                      }
    current_cmd_dict = cmd_label_dict.get(cmd)
    if current_cmd_dict:
        cmd_label = current_cmd_dict["CMD_LABEL"]
        desc = current_cmd_dict["DESC"]
    else:
        cmd_label, desc = None, None
    return(cmd_label,desc)



