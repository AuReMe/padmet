# -*- coding: utf-8 -*-
"""
Description:
    Before running orthology based reconstruction it is necessary to check if the metabolic network 
    and the proteom of the model organism use the same ids for genes (or at least more than a given cutoff). 
    To only check this. Use the 2nd usage.

    If the genes ids are not the same, it is necessary to use a dictionnary of genes ids associating
    the genes ids from the proteom to the genes ids from the metabolic network.

    To create the correct proteom from the dictionnnary, use the 3nd usage
    Finnaly by using the 1st usage, it is possible to:

        1/ Check model_faa and model_metabolic for a given cutoff

        2/ if under the cutoff, convert model_faa to the correct one with dict_ids_file

        3/ if still under, SystemExit()

"""
import re
import itertools
from Bio import SeqIO
import libsbml
from padmet.utils import sbmlPlugin as sp

def check_orthology_input(model_metabolic, model_faa, dict_ids_file, output, verbose, cutoff):
    """
    #TODO
    """
    if model_metabolic is not None:
        if verbose: print("check genes ids model_metablic vs model_faa")
        #if true: more than cutoff% match
        if check_ids(model_metabolic, model_faa, cutoff, verbose):
            return True
        elif dict_ids_file is not None:
            if verbose: print("creating a valid FAA with the dictionnary")
            get_valid_faa(model_faa, dict_ids_file, output)
            if verbose: print("check genes ids")
            if check_ids(model_metabolic, output, cutoff, verbose):
                return True
            else:
                raise SystemExit("Also with the dictionnary, genes ids in the metabolic network and the FAA are not the same")
        else:
                raise SystemExit("Change the cutoff or use an other dictionnary")
                                
    elif dict_ids_file is not None:
        if verbose: print("creating a valid FAA with the dictionnary")
        get_valid_faa(model_faa, dict_ids_file, output)


def check_ids(model_metabolic, model_faa, cutoff, verbose=False):
    """
    check if genes ids of model_metabolic = model_faa for a given cutoff
    faa genes ids are in the first line of each sequence: >GENE_ID ....
    metabolic netowkrs genes ids are in note section, GENE_ASSOCIATION: gene_id-1 or gene_id-2

    Parameters
    ----------
    model_metabolic: str
        path to sbml file
    model_faa: str
        path to fasta faa file
    cutoff: int
        cutoff genes ids from model found in faa
    verbose: bool
        verbose
        
    Returns
    -------
    bool    
        True if same ids, if verbose, print % of genes under cutoff
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_metabolic)
    model = document.getModel()
    document.getNumErrors()
    listOfReactions = model.getListOfReactions()
    #convert to set
    model_metabolic_ids = set(itertools.chain.from_iterable([sp.parseGeneAssoc(geneAssoc) 
    for geneAssoc in (sp.parseNotes(r).get("GENE_ASSOCIATION",[None])[0] for r in listOfReactions)
    if geneAssoc is not None]))
    
    with open(model_faa, "r") as f:
        model_faa_ids = set([record.id for record in SeqIO.parse(f, "fasta")])

    diff_genes = model_metabolic_ids.difference(model_faa_ids)
    try:
        diff_genes_ratio = float(len(diff_genes))/float(len(model_metabolic_ids))
    except ZeroDivisionError:
        raise SystemExit("No genes found in model metabolic")
    #if all model_metabolic_ids are in model_faa_ids
    if diff_genes_ratio == 0:
        if verbose: print("all genes of the model_metabolic are in the model_faa")
        return True
    #if not check if the nb is sup-equal to the cutoff
    elif diff_genes_ratio <= float(1-cutoff):
        if verbose: print("Only %.2f%% genes of the model_metabolic are not in the model_faa" % (diff_genes_ratio*100))
        return True
    else:
        if verbose: 
            print("%s%% genes of the model_metabolic are not in the model_faa" % (diff_genes_ratio*100))
            print(";".join(diff_genes))
        return False

def get_valid_faa(model_faa, dict_ids_file, output):
    """
    create a new faa from the model_faa by converting the gene id with the dict_ids
    dict_ids: line = origin_id new_gene_id, sep = \t

    Parameters
    ----------
    model_faa: str
        path to faa file
    dict_ids_file: str
        path to file containing link old to new ids
    output: str
        path to new faa file
    """
    regex_origin_id = re.compile("^>(\w*)")
    with open(dict_ids_file, 'r') as f:
        dict_ids = dict([(line.split("\t")) for line in f.read().splitlines()])

    with open(output, 'w') as o:
        with open(model_faa, 'r') as f:
            for line in f.readlines():
                line = line.replace('"','')
                if line.startswith(">"):
                    origin_id = regex_origin_id.search(line).group(1)
                    new_gene_id = dict_ids.get(origin_id,origin_id)
                    line = line.replace(origin_id, new_gene_id)
                o.write(line)


