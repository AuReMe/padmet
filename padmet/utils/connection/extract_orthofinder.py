# -*- coding: utf-8 -*-
"""
Description:
    After running orthofinder on n fasta file, read the output file 'Orthogroups.tsv'
    
    Require a folder 'orthology_based_folder' with this archi:
    
        |-- model_a
            -- model_a.sbml
        |-- model_b
            --model_b.sbml

    And the name of the studied organism 'study_id'

    1. Read the orthogroups file, extract orthogroups in dict 'all_orthogroups', and all org names

    2. In orthology folder search for sbml files 'extension = .sbml'

    3. For each models regroup all information in a dict dict_data:
        
        {'study_id': study_id,
        'model_id' : model_id,
        'sbml_template': path to sbml of model',
        'output': path to the output sbml,
        'verbose': bool, if true print information
        }

        The output is by default:
            \output_orthofinder_from_'model_id'.sbml

    4. Store all previous dict_data in a list all_dict_data

    5. iter on dict from all_dict_data and use function dict_data_to_sbml

    Use a dict of data dict_data and dict of orthogroups dict_orthogroup to create sbml files.
    
    dict_data and dict_orthogroup are obtained with fun orthofinder_to_sbml
    
    6./ Read dict_orthogroups and check if model associated to dict_data and study org share orthologue
    
    7./ Read sbml of model, parse all reactions and get genes associated to reaction.
    
    8./ For each reactions:
        
        Parse genes associated to sub part (ex: (gene-a and gene-b) or gene-c) = [(gene-a,gene-b), gene-c]
        
        Check if study org have orthologue with at least one sub part (gene-a, gene-b) or gene-c
        
        if yes: add the reaction to the new sbml and change genes ids by study org genes ids
        
        Create the new sbml file.

::

    usage:
        padmet extract_orthofinder --sbml=FILE/DIR --orthologues=DIR --study_id=STR --output=DIR [--workflow=STR] [-v]
        padmet extract_orthofinder --sbml=DIR --orthogroups=FILE --study_id=STR --output=DIR [--workflow=STR] [-v]

    option:
        -h --help    Show help.
        --sbml=DIR   Folder with sub folder named as models name within sbml file name as model_name.sbml
        --orthogroups=FILE   Output file of Orthofinder run Orthogroups.tsv
        --orthologues=DIR   Output directory of Orthofinder run Orthologues
        --study_id=ID   name of the studied organism
        --workflow=ID   worklow id in ['aureme','aucome']. specific run architecture where to search sbml files
       --output=DIR   folder where to create all sbml output files
        -v   print info
"""
import docopt
import re
import csv
import libsbml
import os

from padmet.utils import sbmlPlugin as sp
from padmet.utils import gbr


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def extract_orthofinder_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    verbose = args["-v"]
    sbml = args["--sbml"]
    orthogroups_file = args["--orthogroups"]
    orthologue_folder = args["--orthologues"]
    output_folder = args["--output"]
    study_id = args["--study_id"]
    workflow = args["--workflow"]
    all_model_sbml = get_sbml_files(sbml, workflow, verbose)
    if orthogroups_file:
        orthogroups_to_sbml(orthogroups_file, all_model_sbml, output_folder, study_id, verbose)
    elif orthologue_folder:
        orthologue_to_sbml(orthologue_folder, all_model_sbml, output_folder, study_id, verbose)


def get_sbml_files(sbml, workflow = None, verbose = False):
    """
    #TODO
    """
    if workflow: workflow = workflow.lower()
    if workflow not in [None, 'aureme','aucome']:
        raise ValueError('--worfklow must be in ["aureme","aucome"]')

    all_model_sbml = {}
    if workflow == 'aureme':
        orthology_based_folder = os.path.join(sbml, "orthology_based_reconstruction")
        if not os.path.exists(orthology_based_folder):
            raise ValueError("folder orthology_based_reconstruction not found in %s, should exist if workflow used is %s" %(sbml, workflow))
        for path, _, files in os.walk(orthology_based_folder):
            fold = os.path.basename(path)
            if fold != "orthofinder_wd":
                for _file in files:
                    if _file == fold+".sbml":
                        all_model_sbml[fold] = os.path.join(path, _file)
    elif workflow == 'aucome':
        annotation_based_sbml_folder = os.path.join(sbml, "annotation_based/SBMLs")
        if not os.path.exists(annotation_based_sbml_folder):
            raise ValueError("folder annotation_based/SBMLs not found in %s, should exist if workflow used is %s" %(sbml, workflow))
        for _file in next(os.walk(annotation_based_sbml_folder))[2]:
            org_id = re.sub("output_pathwaytools_|.sbml", "", _file)
            all_model_sbml[org_id] = os.path.join(annotation_based_sbml_folder, _file)
        model_organisms_folder = os.path.join(sbml, "model_organisms")
        for path, _, files in os.walk(model_organisms_folder):
            org_id = os.path.basename(path)
            for _file in files:
                if _file == org_id+".sbml":
                    all_model_sbml[org_id] = os.path.join(path, _file)
    elif not os.path.isdir(sbml):
        org_id = os.path.splitext(os.path.basename(sbml))[0]
        all_model_sbml[org_id] = sbml
    else:
        for _file in next(os.walk(sbml))[2]:
            org_id = os.path.splitext(_file)[0]
            all_model_sbml[org_id] = os.path.join(sbml, _file)
        
    return all_model_sbml
        
    
    
def orthogroups_to_sbml(orthogroups_file, all_model_sbml, output_folder, study_id, verbose = False):
    """
    After running orthofinder on n fasta file, read the output file 'Orthogroups.tsv'
    Require a folder 'orthology_based_folder' with this archi:
    \model_a
        model_a.sbml
    \model_b
        model_b.sbml
    And the name of the studied organism 'study_id'
    1. Read the orthogroups file, extract orthogroups in dict 'all_orthogroups', and all org names
    2. In orthology folder search for sbml files 'extension = .sbml'
    3. For each models regroup all information in a dict dict_data:
        {'study_id': study_id,
        'model_id' : model_id,
        'sbml_template': path to sbml of model',
        'output': path to the output sbml,
        'verbose': bool, if true print information
        }
        The output is by default: output_orthofinder_from_'model_id'.sbml
    4. Store all previous dict_data in a list all_dict_data
    5. iter on dict from all_dict_data and use function dict_data_to_sbml
    This function will create a sbml from each model and conserve only reactions associated to ortholog genes
    For more information read the doc of func dict_data_to_sbml

    Parameters
    ----------
    orthogroups_file: str
        path of Orthofinder output file 'Orthogroups.tsv'
    orthology_based_folder: str
        path of folder with model's sbml
    output: str
        pathname of the output folder of all sbml extracted
    study_id: str
        name of the studied organism
    verbose: bool
        if True print information
    """
    if verbose:
        print("Parsing Orthofinder output %s" %orthogroups_file)
    #k=orthogroup_id, v = {k = name, v = set of genes}
    dict_orthogroups = {}
    all_orgs = set()
    with open(orthogroups_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter = "\t")
        for row in reader:
            orth_id = row['Orthogroup']
            row.pop('Orthogroup')
            new_dict = dict([(name, set([gene_id.split("_isoform")[0] for gene_id in set(genes.split(","))])) for (name, genes) in list(row.items()) if genes])
            dict_orthogroups[orth_id] = new_dict
            all_orgs.update(new_dict.keys())
    #check all sbml

    if verbose:
        print("Start sbml creation...")
    all_dict_data = []
    if not os.path.exists(output_folder):
        if verbose:
            print("\tCreating folder %s" %output_folder)
        os.makedirs(output_folder)
    all_models = all_orgs - set([study_id])
    for model_id in all_models:
        output_name = "output_orthofinder_from_{0}.sbml".format(model_id)
        output = os.path.join(output_folder, output_name)
        sbml_template = all_model_sbml[model_id]
        dict_data = {'study_id': study_id, 'model_id': model_id, 'sbml_template': sbml_template, 'output': output, 'verbose':verbose}
        all_dict_data.append(dict_data)

    for dict_data in all_dict_data:
        dict_data_to_sbml(dict_data, dict_orthogroups=dict_orthogroups)

def orthologue_to_sbml(orthologue_folder, all_model_sbml, output_folder, study_id, verbose=False):
    """
    After running orthofinder on n fasta file, read the output files in 'Orthologues'
    Require a folder 'orthology_based_folder' with this archi:
    \model_a
        model_a.sbml
    \model_b
        model_b.sbml
    And the name of the studied organism 'study_id'
    1. Read the orthogroups file, extract orthogroups in dict 'all_orthogroups', and all org names
    2. In orthology folder search for sbml files 'extension = .sbml'
    3. For each models regroup all information in a dict dict_data:
        {'study_id': study_id,
        'model_id' : model_id,
        'sbml_template': path to sbml of model',
        'output': path to the output sbml,
        'verbose': bool, if true print information
        }
        The output is by default: output_orthofinder_from_'model_id'.sbml
    4. Store all previous dict_data in a list all_dict_data
    5. iter on dict from all_dict_data and use function dict_data_to_sbml
    This function will create a sbml from each model and conserve only reactions associated to ortholog genes
    For more information read the doc of func dict_data_to_sbml

    Parameters
    ----------
    orthologue_folder: str
        path of Orthofinder output folder 'Orthologues'
    orthology_based_folder: str
        path of folder with model's sbml
    output: str
        pathname of the output folder of all sbml extracted
    study_id: str
        name of the studied organism
    verbose: bool
        if True print information
    """
    if verbose:
        print("Parsing Orthofinder output {0} for {1}".format(orthologue_folder, study_id))
    #k=orthogroup_id, v = {k = name, v = set of genes}
    all_orgs = set()
    all_orthologue_files = []
    for _path, _folders, _files in os.walk(orthologue_folder):
        for _folder in _folders:
            all_orgs.add(_folder.replace("Orthologues_",""))
        for _file in _files:
            _filename = os.path.splitext(_file)[0]
            if _filename.endswith(study_id):
                all_orthologue_files.append(os.path.join(_path,_file))
    dict_orthologues = {}
    for org in all_orgs:
        dict_orthologues[org] = dict()
    #dict_orthologue: k = org_id, v = dict: k = gene_id, v = dict: k = org_id, v = set of gene orthologue id
    for orthologue_file in [i for i in all_orthologue_files if "__v__" in i]:
        with open(orthologue_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            orgs = list(reader.fieldnames)
            orgs.remove('Orthogroup')
            org_A, org_B = orgs
            for row in reader:
                gene_ids_A = [gene_id.split("_isoform")[0] for gene_id in row[org_A].split(", ")]
                gene_ids_B = [gene_id.split("_isoform")[0] for gene_id in row[org_B].split(", ")]
                for gene_id_A in gene_ids_A:
                    if gene_id_A not in dict_orthologues[org_A].keys():
                        dict_orthologues[org_A][gene_id_A] = dict()
                    try:
                        dict_orthologues[org_A][gene_id_A][org_B] = set(gene_ids_B)
                    except KeyError:
                        dict_orthologues[org_A][gene_id_A] = {org_B:set(gene_ids_B)}
                """
                for gene_id_B in gene_ids_B:
                    if gene_id_B not in dict_orthologue[org_B].keys():
                        dict_orthologue[org_B][gene_id_B] = dict()
                    try:
                        dict_orthologue[org_B][gene_id_B][org_A] = set(gene_ids_A)
                    except KeyError:
                        dict_orthologue[org_B][gene_id_B] = {org_A:set(gene_ids_A)}
                """

    if verbose:
        print("Start sbml creation for " + study_id)
    all_dict_data = []
    if not os.path.exists(output_folder):
        if verbose:
            print("\tCreating folder %s" %output_folder)
        os.makedirs(output_folder)        
    all_models = all_orgs - set([study_id])
    for model_id in all_models:
        output_name = "output_orthofinder_from_{0}.sbml".format(model_id)
        output = os.path.join(output_folder, output_name)
        sbml_template = all_model_sbml[model_id]
        dict_data = {'study_id': study_id, 'model_id': model_id, 'sbml_template': sbml_template, 'output': output, 'verbose':verbose}
        all_dict_data.append(dict_data)

    for dict_data in all_dict_data:
        dict_data_to_sbml(dict_data, dict_orthologues=dict_orthologues)


def dict_data_to_sbml(dict_data, dict_orthogroups=None, dict_orthologues=None, strict_match=True):
    """
    Use a dict of data dict_data and dict of orthogroups dict_orthogroup to create sbml files.
    dict_data and dict_orthogroup are obtained with fun orthofinder_to_sbml
    1./ Read dict_orthogroups and check if model associated to dict_data and study org share orthologue
    2./ Read sbml of model, parse all reactions and get genes associated to reaction.
    3./ For each reactions:
        Parse genes associated to sub part (ex: (gene-a and gene-b) or gene-c) = [(gene-a,gene-b), gene-c]
        Check if study org have orthologue with at least one sub part (gene-a, gene-b) or gene-c
        if yes: add the reaction to the new sbml and change genes ids by study org genes ids
    4./ Create the new sbml file.
    
    Parameters
    ----------
    dict_data: dict
        {'study_id': study_id,
        'model_id' : model_id,
        'sbml_template': path to sbml of model',
        'output': path to the output sbml,
        'verbose': bool, if true print information
        }
    dict_orthogroup: dict
        k=orthogroup_id, v = {k = name, v = set of genes}
    verbose: bool
        if True print information
    """
    #dict_data = {'study_name':'', 'o_compare_name': '', sbml_template':'', 'output':''}
    study_id = dict_data['study_id']
    model_id = dict_data['model_id']
    sbml_template = dict_data['sbml_template']
    output = dict_data['output']
    verbose = dict_data.get('verbose')

    if dict_orthogroups:
        if verbose:
            print("*Extracting orthogroups data to create sbml of {0} from {1}".format(study_id, model_id))
    
        #k = gene_id from to_compare, v = list of genes id of study
        sub_dict_orth = {}
        for k in dict_orthogroups.values():
            try:
                all_to_compare_genes = k[model_id]
                all_study_genes = k[study_id]
                for to_compare_gene in all_to_compare_genes:
                    try:
                        sub_dict_orth[to_compare_gene].update(all_study_genes)
                    except KeyError:
                        sub_dict_orth[to_compare_gene] = set(all_study_genes)
            except KeyError:
                pass
                        
        if not sub_dict_orth:
            if verbose:
                print("\t{0} and {1} don't share any ortholgue".format(study_id, model_id))
            return
    elif dict_orthologues:
        if verbose:
            print("*Extracting orthologues data to create sbml of {0} from {1}".format(study_id, model_id))
    
        #k = gene_id from to_compare, v = list of genes id of study
        sub_dict_orth = {}
        for gene_id, gene_dict in dict_orthologues[model_id].items():
            try:
                sub_dict_orth[gene_id] = gene_dict[study_id]
            except KeyError:
                pass
        if not sub_dict_orth:
            if verbose:
                print("\t{0} and {1} don't share any ortholgue".format(study_id, model_id))
            return
    else:
        ValueError("Must give one dict of orthogroups or orthologue")

    reader = libsbml.SBMLReader()
    document_to_compare = reader.readSBML(sbml_template)
    for i in range(document_to_compare.getNumErrors()):
        print(document_to_compare.getError(i).getMessage())
    model_to_compare = document_to_compare.getModel()
    listOfReactions_with_genes = [rxn for rxn in model_to_compare.getListOfReactions()
                                  if sp.parseNotes(rxn).get("GENE_ASSOCIATION",[None])[0]]
    if verbose:
        print("\tSbml of {0} contains {1}/{2} reactions with genes assocation".format(model_id, len(listOfReactions_with_genes), len(model_to_compare.getListOfReactions())))
    dict_rxn_ga = {}
    for rxn in listOfReactions_with_genes:
        ga = sp.parseNotes(rxn)['GENE_ASSOCIATION'][0]
        ga_for_gbr = re.sub(r" or " , "|", ga)
        ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
        ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
        if re.findall("\||&", ga_for_gbr):
            to_compare_ga_subsets = list(gbr.compile_input(ga_for_gbr))
        else:
            ga_for_gbr = re.sub(r"\(|\)" , "", ga_for_gbr)
            to_compare_ga_subsets = [[ga_for_gbr]]
        
        study_ga_subsets = []
        """
        to_compare_ga_subsets = [('a','c','d'),('c',)]
        sub_dict_orth = {'a':['a_a'],'c':['c_c'], 'd':['d_d']}
        """
        for to_compare_subset in to_compare_ga_subsets:
            study_subset = set()
            for gene in to_compare_subset:
                if gene in list(sub_dict_orth.keys()):
                    study_subset.update(sub_dict_orth[gene])
                else:
                    study_subset = set()
                    break
            if study_subset:
                """
                if verbose:
                    print("\t\t{0} == {1}".format(tuple(to_compare_subset), tuple(study_subset)))
                """
                study_ga_subsets.append(study_subset)
        if study_ga_subsets:
            study_ga = " or ".join(["("+" and ".join(subset)+")" for subset in study_ga_subsets])
            if verbose:
                print("\t\tAdding %s" %rxn.id)
                print("\t\tGENE_ASSOCIATION: %s" %(study_ga))
            dict_rxn_ga[rxn.id] = study_ga
    if not dict_rxn_ga:
        if verbose:
            print("\tNo reaction added from {0} to {1} because of missing orthologues".format(model_id, study_id))
        return
    rxn_id_to_remove = set([rxn.id for rxn in model_to_compare.getListOfReactions()]).difference(list(dict_rxn_ga.keys()))
    if verbose:
        print("\tRemoving %s unused reactions" %len(rxn_id_to_remove))
    [model_to_compare.removeReaction(rxn_id) for rxn_id in rxn_id_to_remove]
    cpd_id_to_preserve = set()
    for rxn_id, study_ga in list(dict_rxn_ga.items()):
        rxn = model_to_compare.getElementBySId(rxn_id)
        #update notes
        notes_in_dict = sp.parseNotes(rxn)
        notes_in_dict["GENE_ASSOCIATION"] = [study_ga]
        notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
        for k,v_list in list(notes_in_dict.items()):
            for v in v_list:
                notes += "<p>"+k+": "+v+"</p>"
        notes += "</body>"
        rxn.setNotes(notes)
        cpd_in_rxn = set([p.getSpecies() for p in rxn.getListOfProducts()]).union(\
                         set([r.getSpecies() for r in rxn.getListOfReactants()]))
        cpd_id_to_preserve.update(cpd_in_rxn)
    all_species = [cpd.id for cpd in model_to_compare.getListOfSpecies()]
    [model_to_compare.removeSpecies(cpd_id ) for cpd_id in all_species if cpd_id not in cpd_id_to_preserve]
    new_id = os.path.basename(os.path.splitext(output)[0])    
    model_to_compare.setId(new_id)
    libsbml.writeSBMLToFile(document_to_compare, output)   


