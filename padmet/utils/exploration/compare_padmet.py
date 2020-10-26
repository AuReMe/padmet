# -*- coding: utf-8 -*-
"""
Description:
    #Compare 1-n padmet and create a folder output with files:
    genes.tsv:
        fieldnames = [gene, padmet_a, padmet_b, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [gene-a, 1 (if in padmet_a), 1 (if in padmet_b), rxn-1;rxn-2 (names of reactions associated to gene-a in padmet_a), rxn-2]
    reactions.tsv:
        fieldnames = [reaction, padmet_a, padmet_b, padmet_a_genes_assoc, padmet_b_genes_assoc, padmet_a_formula, padmet_b_formula]
        line = [rxn-1, 1 (if in padmet_a), 1 (if in padmet_b), 'gene-a;gene-b; gene-a, 'cpd-1 + cpd-2 => cpd-3', 'cpd-1 + cpd-2 => cpd-3']
    pathways.tsv:
        fieldnames = [pathway, padmet_a_completion_rate, padmet_b_completion_rate, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [pwy-a, 0.80, 0.30, rxn-a;rxn-b; rxn-a]
    compounds.tsv:
        fieldnames = ['metabolite', padmet_a_rxn_consume, padmet_a_rxn_produce, padmet_b_rxn_consume, padmet_rxn_produce]
        line = [cpd-1, rxn-1,'',rxn-1,'']

::

    usage:
        padmet compare_padmet --padmet=FILES/DIR --output=DIR [--padmetRef=FILE] [--cpu INT] [-v]

    option:
        -h --help    Show help.
        --padmet=FILES/DIR    pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
        --output=DIR    pathname of the output folder
        --padmetRef=FILE    pathanme of the database ref in padmet
        --cpu INT    number of CPU to use in multiprocessing
"""
import docopt
import csv
import os

from multiprocessing import Pool
from padmet.classes import PadmetRef, PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def compare_padmet_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    output = args["--output"]
    verbose = args["-v"]
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    padmet_path = args["--padmet"]
    number_cpu = args["--cpu"]
    compare_padmet(padmet_path, output, padmetRef, verbose, number_cpu)


def extract_information_padmet(file_path, padmetRef, verbose):
    count_elem, dict_genes, dict_rxns, dict_pwys, dict_cpds = {}, {}, {}, {}, {}
    #cout_elem: k = basename_file, v = dict(k' in  [gene,rxn,pwy,metabolite], v = count element)
    #ex: count_elemn = {'file_1':{'genes':1500,'rxn':1000,pathways:1000,'metabolites':1000}}
    #dict_genes: k = gene_id, v = dict(k' = basename_file, v = set of rxn assoc)
    #ex: dict_genes = {'gene_a': {'file_1': set([rxn_1,rxn_2])}}
    #dict_rxns: k = rxn_id, v = dict(k' = basename_file, v' = dict(k'' in ['genes_assoc','formula'])
    #ex: dict_rxns = {'rxn_a': {'file_1': {'genes_assoc': set([gene_a,gene_b]), 'formula':'1x + 2y => 3z', 'source_categories':'', source_tools, sources:''}}
    #dict_pwys: k = pwy_id, , v = dict(k' = basename_file, v' = dict(k'' in ['ratio','rxn_assoc'])
    #ex: dict_pwys = {'pwy_a': {'file_1': {'ratio': 2/9, 'rxn_assoc':'rxn_a;rxn_b'}}}
    #dict_cpds
    padmet = PadmetSpec(file_path)
    basename_file = os.path.basename(file_path).replace(".padmet","")
    log_padmet = ""
    log_padmet += "reading %s\n" %basename_file
    count_elem[basename_file] = {"genes":0,"reactions":0,"pathways":0,"metabolites":0}

    #genes
    all_genes = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "gene"])
    count_elem[basename_file]["genes"] = len(all_genes)
    log_padmet += "\t%s genes...\n" %len(all_genes)
    for gene_id in all_genes:
        #get for each gene, the reactiosn associated to and add in dict_gene: dict_gene[gene_id] = {padmet_file_name: set of rxn assoc}
        rxn_assoc = set([rlt.id_in for rlt in padmet.dicOfRelationOut.get(gene_id,[]) if rlt.type == "is_linked_to"])
        try:
            dict_genes[gene_id][basename_file] = ";".join(rxn_assoc)
        except KeyError:
            dict_genes[gene_id] = {basename_file: ";".join(rxn_assoc)}

    #reactions
    all_rxns = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"])
    count_elem[basename_file]["reactions"] = len(all_rxns)
    log_padmet += "\t%s reactions...\n" %len(all_rxns)
    for rxn_id in all_rxns:
        try:
            dict_rxns[rxn_id][basename_file] = {"genes_associated":"", "formula":""}
        except KeyError:
            dict_rxns[rxn_id] = {basename_file: {"genes_associated":"", "formula":""}}

        #get for each rxn, genes associated to and add in dict_rxn: dict_rxns[rxn_id] = {padmet_file_name: {'genes_associations':set of rxn assoc}}
        genes_assoc = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rxn_id,[]) if rlt.type == "is_linked_to"])
        dict_rxns[rxn_id][basename_file]["genes_associated"] = ";".join(genes_assoc)
        #create formula
        direction = padmet.dicOfNode[rxn_id].misc["DIRECTION"][0]
        if direction == "UNKNOWN":
            direction = " =>/<=> "
        elif direction == "REVERSIBLE":
            direction = " <=> "
        elif direction == "LEFT-TO-RIGHT":
            direction = " => "
        #reactants and products: list of str: ['stoichio [[cpd_id]][cpd_compart]',...]
        reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out+"["+rlt.misc["COMPARTMENT"][0]+"]"
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "consumes"]
        products = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out+"["+rlt.misc["COMPARTMENT"][0]+"]"
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "produces"]
        #join each list with '+' and direction in center
        formula = " + ".join(reactants)+direction+" + ".join(products)
        dict_rxns[rxn_id][basename_file]["formula"] = formula

    #pathways
    all_pwys = set()
    for rxn_id in all_rxns:
        pwy_in = set([rlt.id_out for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"])
        all_pwys.update(pwy_in)
    count_elem[basename_file]["pathways"] = len(all_pwys)
    log_padmet += "\t%s pathways...\n" %len(all_pwys)
    for pwy_id in all_pwys:
        #get for each pwy, the reactions associated to and add in dict_gene: dict_gene[gene_id] = {padmet_file_name: set of rxn assoc}
        try:
            dict_pwys[pwy_id][basename_file] = {"ratio":"", "rxn_associated":""}
        except KeyError:
            dict_pwys[pwy_id] = {basename_file: {"ratio":"", "rxn_associated":""}}
        in_rxns = set([rlt.id_in for rlt in padmet.dicOfRelationOut.get(pwy_id,[]) if rlt.type == "is_in_pathway"])
        if padmetRef:
            all_rxns = set([rlt.id_in for rlt in padmetRef.dicOfRelationOut.get(pwy_id,[]) if rlt.type == "is_in_pathway"])
            dict_pwys[pwy_id][basename_file]["ratio"] = str(len(in_rxns))+"/"+str(len(all_rxns))
        else:
            dict_pwys[pwy_id][basename_file]["ratio"] = "NA"
        dict_pwys[pwy_id][basename_file]["rxn_associated"] = ";".join(in_rxns)

    #metabolites
    all_cpd = set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type in ["consumes","produces"]])
    count_elem[basename_file]["metabolites"] = len(all_cpd)
    log_padmet += "\t%s metabolites...\n" %len(all_cpd)
    for cpd_id in all_cpd:
        #get for each cpd, the reactions consuming/producing the cpd
        rxn_consume = set([rlt.id_in for rlt in padmet.dicOfRelationOut[cpd_id] if rlt.type == "consumes"])
        if rxn_consume:
            rxn_consume = ";".join(rxn_consume)
        else:
            rxn_consume = ""
        rxn_produce = set([rlt.id_in for rlt in padmet.dicOfRelationOut[cpd_id] if rlt.type == "produces"])
        if rxn_produce:
            rxn_produce = ";".join(rxn_produce)
        else:
            rxn_produce = ""
        try:
            dict_cpds[cpd_id][basename_file] = {"rxn_consume":rxn_consume, "rxn_produce":rxn_produce}
        except KeyError:
            dict_cpds[cpd_id] = {basename_file: {"rxn_consume":rxn_consume, "rxn_produce":rxn_produce}}

    if verbose:
        print(log_padmet)
    return dict_genes, dict_rxns, dict_pwys, dict_cpds


def merge_dicts(element_dict, tmp_dict):
    for gene_id in tmp_dict:
        if gene_id in element_dict:
            for basename_file in tmp_dict[gene_id]:
                element_dict[gene_id][basename_file] = tmp_dict[gene_id][basename_file]
        else:
            element_dict[gene_id] = {}
            for basename_file in tmp_dict[gene_id]:
                element_dict[gene_id][basename_file] = tmp_dict[gene_id][basename_file]

    return element_dict


def compare_padmet(padmet_path, output, padmetRef = None, verbose = False, number_cpu = None):
    """
    #Compare 1-n padmet and create a folder output with files:
    genes.tsv:
        fieldnames = [gene, padmet_a, padmet_b, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [gene-a, 1 (if in padmet_a), 1 (if in padmet_b), rxn-1;rxn-2 (names of reactions associated to gene-a in padmet_a), rxn-2]
    reactions.tsv:
        fieldnames = [reaction, padmet_a, padmet_b, padmet_a_genes_assoc, padmet_b_genes_assoc, padmet_a_formula, padmet_b_formula]
        line = [rxn-1, 1 (if in padmet_a), 1 (if in padmet_b), 'gene-a;gene-b; gene-a, 'cpd-1 + cpd-2 => cpd-3', 'cpd-1 + cpd-2 => cpd-3']
    pathways.tsv:
        fieldnames = [pathway, padmet_a_completion_rate, padmet_b_completion_rate, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [pwy-a, 0.80, 0.30, rxn-a;rxn-b; rxn-a]
    compounds.tsv:
        fieldnames = ['metabolite', padmet_a_rxn_consume, padmet_a_rxn_produce, padmet_b_rxn_consume, padmet_rxn_produce]
        line = [cpd-1, rxn-1,'',rxn-1,'']

    Parameters
    ----------
    padmet_path: str
        pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
    output: str
        pathname of the output folder
    padmetRef: padmet.classes.PadmetRef
        padmet containing the database of reference, need to calculat pathway completion rate
    verbose: bool
        if True print information
    
    """
    dict_genes, dict_rxns, dict_pwys, dict_cpds = {}, {}, {}, {}
    if not os.path.exists(output):
        if verbose: print("Creating %s" %output)
        os.makedirs(output)
    else:
        if verbose: print("%s already exist, old comparison output folders will be overwritten" %output)

    if os.path.isdir(padmet_path):
        all_files = [os.path.join(padmet_path, f) for f in next(os.walk(padmet_path))[2]]
    else:
        all_files = padmet_path.split(",")

    if len(all_files) < 2:
        raise ValueError("You must specify at least 2 files in order to make a comparison")
    if verbose:
        print(("%s padmet files to compare:" %len(all_files)))
        for f in all_files:
            print("\t%s" %os.path.basename(f))

    if number_cpu:
        try:
            number_cpu_to_use = int(number_cpu)
        except ValueError:
            raise ValueError('The number of CPU must be an integer.')
    else:
        number_cpu_to_use = 1
    compare_pool = Pool(processes=number_cpu_to_use)

    multiprocessin_datas = []
    for padmet_file in all_files:
        multiprocessin_datas.append([padmet_file, padmetRef, verbose])
    result_rxns_dicts = compare_pool.starmap(extract_information_padmet, multiprocessin_datas)

    for result_rxns_dict in result_rxns_dicts:
        tmp_dict_genes = result_rxns_dict[0]
        tmp_dict_rxns = result_rxns_dict[1]
        tmp_dict_pwys = result_rxns_dict[2]
        tmp_dict_cpds = result_rxns_dict[3]
        dict_genes = merge_dicts(dict_genes, tmp_dict_genes)
        dict_rxns = merge_dicts(dict_rxns, tmp_dict_rxns)
        dict_pwys = merge_dicts(dict_pwys, tmp_dict_pwys)
        dict_cpds = merge_dicts(dict_cpds, tmp_dict_cpds)

    compare_pool.close()
    compare_pool.join()

    #create files
    all_basename_files = [os.path.basename(file_path).replace(".padmet","") for file_path in all_files]
    #genes
    #gene file header: gene_id, base_file_1, base_file_n, base_file_1_rxn_assoc (sep=;), base_file_n_rxn_assoc (sep=;)
    genes_file = os.path.join(output,"genes.tsv")
    if verbose: print("creating %s" %genes_file)
    with open(genes_file, 'w') as csvfile:
        fieldnames = ['gene'] + all_basename_files + [i+"_rxn_assoc (sep=;)" for i in all_basename_files]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for gene_id, dic_basename_rxn_assoc in list(dict_genes.items()):
            dict_row = {'gene': gene_id}
            for basename_file, rxn_assoc in list(dic_basename_rxn_assoc.items()):
                dict_row.update({basename_file : 1, basename_file+"_rxn_assoc (sep=;)": rxn_assoc})
            writer.writerow(dict_row)

    #reactions
    #reactions file header: rxn_id, base_file_1, base_file_n, base_file_1_genes_assoc (sep=;), base_file_n_genes_assoc (sep=;), base_file_1_formula, base_file_n_formula
    rxns_file = os.path.join(output,"reactions.tsv")
    if verbose: print("creating %s" %rxns_file)
    with open(rxns_file, 'w') as csvfile:
        fieldnames = ['reaction'] + all_basename_files + [i+"_genes_assoc (sep=;)" for i in all_basename_files] + [i+"_formula" for i in all_basename_files]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for rxn_id, dict_basename_data in list(dict_rxns.items()):
            dict_row = {'reaction': rxn_id}
            for basename_file, rxn_data in list(dict_basename_data.items()):
                dict_row.update({basename_file: 1, basename_file+"_genes_assoc (sep=;)": rxn_data["genes_associated"], basename_file+"_formula": rxn_data["formula"]})
            for basename_file in all_basename_files:
                if basename_file not in dict_row:
                    dict_row.update({basename_file: 0})
            writer.writerow(dict_row)

    #pathways
    #pathways file header: pwy, base_file_1_rate, base_file_n_rate, base_file_1_rxn_assoc (sep=;), base_file_n_rxn_assoc (sep=;)
    pwys_file = os.path.join(output,"pathways.tsv")
    if verbose: print("creating %s" %pwys_file)
    with open(pwys_file, 'w') as csvfile:
        fieldnames = ['pathway'] + [i+"_completion_rate" for i in all_basename_files] + [i+"_rxn_assoc (sep=;)" for i in all_basename_files]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for pwy_id, dict_basename_data in list(dict_pwys.items()):
            dict_row = {'pathway': pwy_id}
            for basename_file, pwy_data in list(dict_basename_data.items()):
                dict_row.update({basename_file+"_completion_rate": pwy_data["ratio"], basename_file+"_rxn_assoc (sep=;)": pwy_data["rxn_associated"]})
            writer.writerow(dict_row)

    #metabolites
    #metabolites file header: cpd, base_file_1, base_file_n
    cpds_file = os.path.join(output,"metabolites.tsv")
    if verbose: print(("creating %s" %cpds_file))
    with open(cpds_file, 'w') as csvfile:
        fieldnames = ['metabolite'] + [i+"_rxn_consume" for i in all_basename_files] + [i+"_rxn_produce" for i in all_basename_files]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for cpd_id, dict_basename_data in dict_cpds.items():
            dict_row = {'metabolite': cpd_id}
            for basename_file, cpd_data in dict_basename_data.items():
                dict_row.update({basename_file+"_rxn_consume": cpd_data["rxn_consume"], basename_file+"_rxn_produce": cpd_data["rxn_produce"]})
            writer.writerow(dict_row)


