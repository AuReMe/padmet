# -*- coding: utf-8 -*-
"""
Description:

    Read a PGDB folder (from BIOCYC/PATHWAYTOOLS) and create a padmet.
    1./ To create a padmet without any genes information extracted use the first usage with:
        pgdb: path to pgdb folder
        output: path to the padmet to create
        version: to specify the version of the pgdb (20.0, 22.0)
        db: to sepcify the name of the database (METACYC, ECOCYC, ...)
        enhance: to also read the file metabolic-reaction.xml and add the to the padmet
    2./ To create a padmet and add only reactions from pgdb if they are in padmetRef specifie.
        Copy information of the reaction not from the pgdb but from the padmetRef.
        This allow to uniform reaction to the same version of metacyc represented in the padmetRef
        For example, in some case 2 pgdb from different version can contain different information for a same reaction,pathway...
        In this case use:
            padmetRef: path to the padmet of reference
    3./ To create a padmet wth genes information extracted use:
        extract-gene
    3.1/ To remove from the final padmet all reactions without genes associated use:
        no-orphan
    4./ To read the metabolic-reaction.xml file, a sbml with some missing reactions in PGDB use:
        enhance

    For more information of the parsing process read information below.

    
    
    classes.dat:
    For each class:
    create new node / class = class
    UNIQUE-ID (1) => node.id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name synonyms)
     
    compounds.dat:
    for each compound:
    create new node / class = compound
    UNIQUE-ID (1) => node.id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
    MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
    SMILES (0-1) => node.misc()['SMILES'] = SMILES
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
    DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)
    
    proteins.dat:
    for each protein:
    create new node / class = protein
    UNIQUE-ID (1) => node.id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
    MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
    SMILES (0-1) => node.misc()['SMILES'] = SMILES
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
    DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)
    SPECIES (0-1) => for each, check or create new node class, create rlt (node is_in_species class)
    
    reactions.dat:
    for each reaction:
    create new node / class = reaction + node.misc()["DIRECTION"] = "UNKNOWN" by default
    UNIQUE-ID (1) => node.id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    EC-NUMBER (0-n) => node.Misc['EC-NUMBER'] = EC-NUMBER
    REACTION-DIRECTION (0-1) => node.Misc['DIRECTION'] = reaction-direction, if REVERSIBLE, else: LEFT-TO-RIGHT
    RXN-LOCATIONS (0,n) => node.misc['COMPARTMENT'] = rxn-location
    TYPES (0-n) => check or create new node class, create rlt (node.id is_a_class types's_node.id)
    DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref's_node.id)
    SYNONYMS (0-n) => create new node name, create rlt (node has_name name's_node.id)
    --
    for LEFT and RIGHT, also check 2 next lines if info about 'coefficient' or 'compartment'
    defaut value: coefficient/stoichiometry = 1, compartment = unknown
    also check if the direction is 'RIGHT-TO-LEFT', if yes, inverse consumes and produces relations
    then change direction to 'LEFT-TO-RIGHT'
    LEFT (1-n) => create rlt (node.id consumes left's_node.id)
    RIGHT (1-n) => create rlt (node.id produces right's_node.id)
    
    enzrxns.dat:
    for each association enzyme/reaction:
    create new rlt / type = catalyses
    ENZYME (1) => stock enzyme as 'enzyme catalyses'
    REACTION (1-n) => for each reaction after, create relation 'enzyme catalyses reaction'
    
    pathways.dat:
    for each pathway:
    create new node / class = pathway
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    TYPES (0-n) => check or create new node class, create rlt (node is_a_class types)
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref)
    SYNONYMS (0-n) => create new node name, create rlt (node has_name name)
    IN-PATHWAY (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)
    REACTION-LIST (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)

::

    usage:
        padmet pgdb_to_padmet --pgdb=DIR --output=FILE [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]
        padmet pgdb_to_padmet --pgdb=DIR --output=FILE --extract-gene [--no-orphan] [--keep-self-rxn] [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]

    options:
        -h --help     Show help.
        --version=V    Xcyc version [default: N.A].
        --db=ID    Biocyc database corresponding to the pgdb (metacyc, ecocyc, ...) [default: N.A].
        --output=FILE    padmet file corresponding to the DB.
        --pgdb=DIR    directory containg all the .dat files of metacyc (data).
        --padmetRef=FILE    padmet of reference.
        --source=STR    Tag associated to the source of the reactions, used to ensure traceability [default: GENOME].
        --enhance    use the metabolic-reactions.xml file to enhance the database.
        --extract-gene    extract genes from genes_file (use if its a specie's pgdb, if metacyc, do not use).
        --no-orhpan    remove reactions without gene associaiton (use if its a specie's pgdb, if metacyc, do not use).
        --keep-self-rxn    remove reactions with no reactants (use if its a specie's pgdb, if metacyc, do not use).
        -v   print info.
"""
import docopt
import libsbml
import re
import os

from padmet.classes import PadmetRef, Node, Relation, instantiate_padmet
import padmet.utils.sbmlPlugin as sbmlPlugin


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def pgdb_to_padmet_cli(command_args):
    #parsing args
    args = docopt.docopt(__doc__, argv=command_args)
    version = args["--version"]
    db = args["--db"]
    source = args["--source"]
    output = args["--output"]
    pgdb_folder = args["--pgdb"]
    enhanced_db = args["--enhance"]
    extract_gene = args["--extract-gene"]
    no_orphan = args["--no-orphan"]
    keep_self_producing_rxn = args["--keep-self-rxn"]
    padmetRef_file = args["--padmetRef"]
    verbose = args["-v"]

    if keep_self_producing_rxn:
        no_self_producing_rxn = False
    else:
        no_self_producing_rxn = True

    from_pgdb_to_padmet(pgdb_folder=pgdb_folder, db=db , version=version, source=source, extract_gene=extract_gene,
                                                no_orphan=no_orphan, no_self_producing_rxn=no_self_producing_rxn, enhanced_db=enhanced_db,
                                                padmetRef_file=padmetRef_file, verbose=verbose, output_file=output)


def from_pgdb_to_padmet(pgdb_folder, db='MetaCyc', version='NA', source='GENOME', extract_gene=False, no_orphan=False, no_self_producing_rxn=True, enhanced_db=False, padmetRef_file=None, verbose=False, output_file=None):
    """
    Parameters
    ----------
    pgdb_folder: str
        path to pgdb
    db: str
        pgdb name, default is 'NA'
    version: str
        pgdb version, default is 'NA'
    source: str
        tag reactions for traceability, default is 'GENOME'
    extract_gene: bool
        if true extract genes information
    no_orphan: bool
        if true, remove reactions without genes associated
    no_self_producing_rxn: bool
        if true, remove reactions with no reactants (auto-producing reactions)
    enhanced_db: bool
        if true, read metabolix-reactions.xml sbml file and add information in final padmet
    padmetRef_file: str
        path to padmetRef corresponding to metacyc in padmet format
    verbose: bool
        if True print information
    output_file: str
        pathname of padmet output file

    Returns
    -------
    padmet.padmetRef:
        padmet instance with pgdb within pgdb data
    """
    global regex_purge, regex_xref, list_of_relation, def_compart_in, def_compart_out
    regex_purge = re.compile(r"<.*?>|\|")
    regex_xref = re.compile(r'^\((?P<DB>\S*)\s*"(?P<ID>\S*)"')
    list_of_relation = []
    def_compart_in = "c"
    def_compart_out = "e"
    #parsing args
    source = source.upper()

    if output_file:
        padmet_id = os.path.splitext(os.path.basename(output_file))[0]
    else:
        padmet_id = None

    if not os.path.exists(pgdb_folder):
        raise FileNotFoundError("No PGDB folder (--pgdb/pgdb_folder) accessible at " + pgdb_folder)

    classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file = \
    [os.path.join(pgdb_folder,_file) for _file in ["classes.dat", "compounds.dat", "proteins.dat", "reactions.dat", "enzrxns.dat", "pathways.dat"]]
    if enhanced_db:
        metabolic_reactions = os.path.join(pgdb_folder,"metabolic-reactions.xml")
    else:
        metabolic_reactions = None
    if extract_gene:
        genes_file = os.path.join(pgdb_folder,"genes.dat")
    else:
        genes_file = None

    if padmetRef_file:
        padmetRef = PadmetRef(padmetRef_file)

        padmet = instantiate_padmet("PadmetSpec", padmetRef_file, padmet_id, db, version, verbose)

        with open(reactions_file, 'r') as f:
            rxns_id = [line.split(" - ")[1] for line in f.read().splitlines() if line.startswith("UNIQUE-ID")]
        count = 0
        for rxn_id in rxns_id:
            count += 1
            if verbose: print("%s/%s Copy %s" %(count, len(rxns_id), rxn_id))
            try:
                padmet.copyNode(padmetRef, rxn_id)
                reconstructionData_id = rxn_id+"_reconstructionData_"+source
                if reconstructionData_id in list(padmet.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(rxn_id, source))
                reconstructionData = {"SOURCE":[source],"TOOL":["PATHWAYTOOLS"],"CATEGORY":["ANNOTATION"]}
                reconstructionData_rlt = Relation(rxn_id,"has_reconstructionData",reconstructionData_id)
                padmet.dicOfNode[reconstructionData_id] = Node("reconstructionData", reconstructionData_id, reconstructionData)
                padmet._addRelation(reconstructionData_rlt)

            except TypeError:
                print("%s not in padmetRef" %(rxn_id))

        if verbose: print("parsing compounds")
        compounds = compounds_parser(compounds_file, padmet, verbose)
        if verbose: print("parsing classes")
        id_classes = classes_parser(classes_file, padmet, verbose)
        if verbose: print("parsing rnas")
        rnas = rnas_parser(os.path.join(pgdb_folder,'rnas.dat'), padmet, verbose)

        if extract_gene:
            if verbose: print("parsing genes")
            map_gene_ids = genes_parser(genes_file, padmet, verbose)
            if verbose: print("parsing proteins")
            dict_protein_gene_id = proteins_parser(proteins_file, padmet, compounds, rnas, id_classes, verbose)
            mapped_dict_protein_gene_id = map_gene_id(dict_protein_gene_id, map_gene_ids)
            if verbose: print("parsing association enzrxns")
            enzrxns_parser(enzrxns_file, padmet, mapped_dict_protein_gene_id, source, verbose)

    else:
        padmet = instantiate_padmet("PadmetRef", None, padmet_id, db, version, verbose)

        if verbose: print("parsing classes")
        id_classes = classes_parser(classes_file, padmet, verbose)
    
        if verbose: print("parsing compounds")
        compounds = compounds_parser(compounds_file, padmet, verbose)

        if verbose: print("parsing rnas")
        rnas = rnas_parser(os.path.join(pgdb_folder,'rnas.dat'), padmet, verbose)
    
        if verbose: print("parsing reactions")
        reactions_parser(reactions_file, padmet, extract_gene, source, verbose)

        if verbose: print("parsing pathways")
        pathways_parser(pathways_file, padmet, verbose)

        if extract_gene:
            if verbose: print("parsing genes")
            map_gene_ids = genes_parser(genes_file, padmet, verbose)
            if verbose: print("parsing proteins")
            dict_protein_gene_id = proteins_parser(proteins_file, padmet, compounds, rnas, id_classes, verbose)
            mapped_dict_protein_gene_id = map_gene_id(dict_protein_gene_id, map_gene_ids)
            if verbose: print("parsing association enzrxns")
            enzrxns_parser(enzrxns_file, padmet, mapped_dict_protein_gene_id, source, verbose)
    
        if metabolic_reactions is not None:
            if verbose: print("enhancing db from metabolic-reactions.xml")
            padmet = enhance_db(metabolic_reactions, padmet, extract_gene, verbose)
    
    for rlt in list_of_relation:
        try:
            padmet.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmet.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmet.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmet.dicOfRelationOut[rlt.id_out] = [rlt]

    if extract_gene and no_orphan:
        all_reactions = [node for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        rxn_to_del = [r for r in all_reactions if not any([rlt for rlt in padmet.dicOfRelationIn[r.id] if rlt.type == "is_linked_to"])]
        for rxn in rxn_to_del:
            if 'SPONTANEOUS' not in rxn.misc:
                padmet.delNode(rxn.id)
        if verbose:
            print("%s/%s orphan reactions (without gene association) deleted" %(len(rxn_to_del), len(all_reactions)))
        all_genes_linked = set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type == "is_linked_to"])
        all_genes = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "gene"])
        count = 0
        for gene_id in [g for g in all_genes if g not in all_genes_linked]:
            count += 1
            #if verbose: print("Removing gene without gene assoc %s" %gene_id)
            padmet.dicOfNode.pop(gene_id)
        if verbose:
            print("%s/%s orphan genes (not linked to any reactions) deleted" %(count, len(all_genes)))
    if no_self_producing_rxn:
        rxns = [node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        for rxn_id in rxns:
            cp_rlts = set([rlt.type for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]])
            if len(cp_rlts) == 1:
                print("rxn only consume or produce, transport ???: %s" %rxn_id)
                padmet.delNode(rxn_id)

    if output_file:
        padmet.generateFile(output_file)

    return padmet

def classes_parser(filePath, padmet, verbose = False):
    """
    from class.dat: get for each class, the UNIQUE-ID, COMMON-NAME, TYPES, SYNONYMS, DBLINKS
    Create a class node with node.id = UNIQUE-ID,  node.misc = {COMMON-NAME:[COMMON-NAMES]}
    - For each types:
    A type is in fact a class. this information is stocked in padmet as: is_a_class relation btw a node and a class_node
    check if the type is already in the padmet
    if not create a new class_node (var: subClass) with subClass_node.id = type
    Create a relation current node is_a_class type
    - For each Synonyms:
    this information is stocked in padmet as: has_name relation btw a node and a name_node
    create a new name_node with name_node.id = class_id+"_names" and name_node.misc = {LABEL:[synonyms]}
    Create a relation current node has_name name_node.id
    - For each DBLINKS:
    DBLINKS is parsed with regex_xref to get the db and the id
    this information is stocked in padmet as: has_xref relation btw a node and a xref_node
    create a new xref_node with xref_node.id = class_id+"_xrefs" and xref_node.misc = {db:[id]}
    Create a relation current node has_xref xref_node.id

    Parameters
    ----------
    filePath: str
        path to classes.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """
    id_classes = []
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    id_classes.append(current_id)
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

    count = 0
    nb_classes = str(len(list(dict_data.keys())))
    for class_id, dict_values in dict_data.items():
        count += 1
        if verbose:
            print("\r%s/%s" %(count, nb_classes), end="", flush=True)
            #print(class_id)
        class_node = Node("class", class_id)
        padmet.dicOfNode[class_id] = class_node
        try:
            class_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types, class_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, class_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, class_id, padmet)
        except KeyError:
            pass
    if verbose: print("")
    return id_classes


def reactions_parser(filePath, padmet, extract_gene, source, verbose = False):
    """
    from reaction.dat: get for each reaction, the UNIQUE-ID, COMMON-NAME, TYPES, SYNONYMS, DBLINKS
    Create a reaction node with node.id = UNIQUE-ID,  node.misc = {COMMON-NAME:[COMMON-NAMES]}
    - For each types:
    A type is in fact a class. this information is stocked in padmet as: is_a_class relation btw a node and a class_node
    check if the type is already in the padmet
    if not create a new class_node (var: subClass) with subClass_node.id = type
    Create a relation current node is_a_class type
    - For each Synonyms:
    this information is stocked in padmet as: has_name relation btw a node and a name_node
    create a new name_node with name_node.id = reaction_id+"_names" name_node.misc = {LABEL:[synonyms]}
    Create a relation current node has_name name_node.id
    - For each DBLINKS:
    DBLINKS is parsed with regex_xref to get the db and the id
    this information is stocked in padmet as: has_xref relation btw a node and a xref_node
    create a new xref_node with xref_node.id = reaction_id+"_xrefs" and xref_node.misc = {db:[id]}
    Create a relation current node has_xref xref_node.id

    Parameters
    ----------
    filePath: str
        path to reactions.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = [line for line in f.read().splitlines() if not line.startswith("#") and not line == "//"]
        index = -1        
        for line in data:
            index += 1
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME", "EC-NUMBER", "REACTION-DIRECTION",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

                elif attrib in ["LEFT", "RIGHT"]:
                    #set default values
                    compartment = def_compart_in
                    stoichiometry = "1"
                    #check if information about stoechiometry and compartment in line + 1 and line + 2
                    try:
                        first_next_line = data[index+1]
                    #last line of the file, just add the LEFT/RIGHT
                    except IndexError:
                        first_next_line = ""
                    try:
                        second_next_line = data[index+2]
                    except IndexError:
                        second_next_line = ""

                    if first_next_line.startswith("^COEFFICIENT"):
                        stoichiometry = first_next_line.split(" - ")[1]
                        if second_next_line.startswith("^COMPARTMENT"):
                            compartment = second_next_line.split(" - ")[1]
                    if first_next_line.startswith("^COMPARTMENT"):
                        compartment = first_next_line.split(" - ")[1]
                        if second_next_line.startswith("^COEFFICIENT"):
                            stoichiometry = second_next_line.split(" - ")[1]                            
                    
                    #delete all tags
                    compartment = re.sub(regex_purge, "", compartment)
                    stoichiometry = re.sub(regex_purge, "", stoichiometry)
                    try:
                        dict_data[current_id][attrib].append((value, stoichiometry, compartment))
                    except KeyError:
                        dict_data[current_id][attrib] = [(value, stoichiometry, compartment)]
                elif 'SPONTANEOUS' in attrib:
                    if value == 'T':
                        dict_data[current_id]['SPONTANEOUS'] = value

    count = 0
    nb_rxn = str(len(list(dict_data.keys())))
    for rxn_id, dict_values in dict_data.items():
        if "LEFT" in list(dict_values.keys()) or "RIGHT" in list(dict_values.keys()):
            count += 1
            if verbose:
                print("\r%s/%s: %s" %(count, nb_rxn, rxn_id), end="", flush=True)
                #print(rxn_id)
            rxn_node = Node("reaction", rxn_id)
            padmet.dicOfNode[rxn_id] = rxn_node
            try:
                rxn_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
            except KeyError:
                pass
            try:
                rxn_node.misc["EC-NUMBER"] = dict_values["EC-NUMBER"]
            except KeyError:
                pass
            try:
                rxn_dir = dict_values["REACTION-DIRECTION"][0]
                if rxn_dir == "REVERSIBLE":
                    rxn_node.misc["DIRECTION"] = ["REVERSIBLE"]
                elif "LEFT-TO-RIGHT" in rxn_dir:
                    #if:LEFT-TO-RIGHT, IRREVERSIBLE-LEFT-TO-RIGHT, PHYSIOL-RIGHT-TO-LEFT
                    rxn_node.misc["DIRECTION"] = ["LEFT-TO-RIGHT"]
                elif "RIGHT-TO-LEFT" in rxn_dir:
                    #Temporarily set direaction as RIGHT-TO-LEFT
                    #then, RIGHT' metabolites will be LEFT and LEFT -> RIGHT
                    #To finish set back DIRECTION to LEFT-TO-RIGHT
                    rxn_node.misc["DIRECTION"] = ["RIGHT-TO-LEFT"] 
            except KeyError:
                rxn_node.misc["DIRECTION"] = ["REVERSIBLE"]
            try:
                rxn_node.misc["SPONTANEOUS"] = dict_values["SPONTANEOUS"]
            except KeyError:
                pass
            """
            try:
                rxn_node.misc["COMPARTMENT"] = dict_values["RXN-LOCATIONS"]
            except KeyError:
                pass
            """
            try:
                types = dict_values["TYPES"]
                _setType(types, rxn_id, padmet)
            except KeyError:
                pass
            try:
                syns = dict_values["SYNONYMS"]
                _setSyns(syns, rxn_id, padmet)
            except KeyError:
                pass
            try:
                xrefs = dict_values["DBLINKS"]
                _setXrefs(xrefs, rxn_id, padmet)
            except KeyError:
                pass
            if extract_gene:
                reconstructionData_id = rxn_id+"_reconstructionData_"+source
                if reconstructionData_id in list(padmet.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(rxn_id, source))
                reconstructionData = {"SOURCE":[source],"TOOL":["PATHWAYTOOLS"],"CATEGORY":["ANNOTATION"]}
                reconstructionData_rlt = Relation(rxn_id,"has_reconstructionData",reconstructionData_id)
                padmet.dicOfNode[reconstructionData_id] = Node("reconstructionData", reconstructionData_id, reconstructionData)
                list_of_relation.append(reconstructionData_rlt)
            try:
                reactants_data = dict_values["LEFT"]
                for reactant_id, stoichiometry, compartment in reactants_data:
                    if compartment == "CCO-OUT":
                        compartment = def_compart_out
                    else:
                        compartment = def_compart_in
                    try:
                        reactant_node = padmet.dicOfNode[reactant_id]
                    except KeyError:
                        reactant_node = Node("compound", reactant_id)
                        padmet.dicOfNode[reactant_id] = reactant_node
    
                    #if the reaction direction was set to RIGHT-TO-LEFT, then this compound is in fact a product
                    if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                        produces_rlt = Relation(rxn_id, "produces", reactant_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(produces_rlt)
                    else:
                        consumes_rlt = Relation(rxn_id, "consumes", reactant_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(consumes_rlt)
            except KeyError:
                pass
            try:
                products_data = dict_values["RIGHT"]
                for product_id, stoichiometry, compartment in products_data:
                    if compartment == "CCO-OUT":
                        compartment = def_compart_out
                    else:
                        compartment = def_compart_in
                    try:
                        product_node = padmet.dicOfNode[product_id]
                    except KeyError:
                        product_node = Node("compound", product_id)
                        padmet.dicOfNode[product_id] = product_node

                    #if the reaction direction was set to RIGHT-TO-LEFT, then this compound is in fact a reactant
                    if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                        consumes_rlt = Relation(rxn_id, "consumes", product_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(consumes_rlt)
                    else:
                        produces_rlt = Relation(rxn_id, "produces", product_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(produces_rlt)
            except KeyError:
                pass
            if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                rxn_node.misc["DIRECTION"] = ["LEFT-TO-RIGHT"]
    if verbose: print("")

def pathways_parser(filePath, padmet, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to pathways.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """        
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                #if len of value is 0 then ValueError raised
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME", "TAXONOMIC-RANGE",\
                "TYPES", "SYNONYMS", "DBLINKS", "IN-PATHWAY",\
                    "REACTION-LIST", "REACTION-LAYOUT"]:
                    if attrib in dict_data[current_id]:
                        dict_data[current_id][attrib].append(value)
                    else:
                        dict_data[current_id][attrib] = [value]

    count = 0
    nb_pathways = str(len(list(dict_data.keys())))
    for pathway_id, dict_values in dict_data.items():
        count += 1
        if verbose:
            print("\r%s/%s" %(count, nb_pathways), end="", flush=True)
            #print(pathway_id)
        pathway_node = Node("pathway", pathway_id)
        padmet.dicOfNode[pathway_id] = pathway_node
        try:
            pathway_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            pathway_node.misc["TAXONOMIC-RANGE"] = dict_values["TAXONOMIC-RANGE"]
        except KeyError:
            pass

        if "REACTION-LAYOUT" in dict_values:
            reaction_layout_regex = r'\((?P<reaction_id>\S*)\s\(\:LEFT-PRIMARIES[\s]*(?P<reaction_left>[^\(]*)\)\s\(\:DIRECTION\s(?P<reaction_direction>\S*)\)\s\(:RIGHT-PRIMARIES[\s]*(?P<reaction_right>[^\(]*)\)\)'
            reactions = dict_values["REACTION-LAYOUT"]
            reactions_compounds = []
            for reaction in reactions:
                regex_xref = re.match(reaction_layout_regex, reaction)
                if regex_xref is None:
                    print(reaction)
                reaction_id = regex_xref.groupdict()['reaction_id']
                reaction_direction = regex_xref.groupdict()['reaction_direction']
                if reaction_direction == ':L2R':
                    reaction_reactant = regex_xref.groupdict()['reaction_left'].split(' ')
                    reaction_product = regex_xref.groupdict()['reaction_right'].split(' ')
                elif reaction_direction == ':R2L':
                    reaction_reactant = regex_xref.groupdict()['reaction_right'].split(' ')
                    reaction_product = regex_xref.groupdict()['reaction_left'].split(' ')

                reactions_compounds.append((reaction_id, reaction_reactant, reaction_product))

            reactions_ordered = {}
            reactions_added = []

            pathway_reactants = [reactant for reaction in reactions_compounds for reactant in reaction[1]]
            pathway_products = [product for reaction in reactions_compounds for product in reaction[2]]

            pathway_input_compounds = list(set(pathway_reactants) - set(pathway_products))
            pathway_output_compounds = list(set(pathway_products) - set(pathway_reactants))

            if pathway_input_compounds != []:
                pathway_node.misc["INPUT-COMPOUNDS"] = [','.join(pathway_input_compounds)]
            if pathway_output_compounds != []:
                pathway_node.misc["OUTPUT-COMPOUNDS"] = [','.join(pathway_output_compounds)]

            def reaction_order(compound, reaction_ordered, reactions_compounds, alreayd_known_compounds):
                for reaction in reactions_compounds:
                    if compound in reaction[1]:
                        reaction_ordered.append(reaction[0])
                        alreayd_known_compounds.append(compound)
                        for reaction_product in reaction[2]:
                            if reaction_product not in alreayd_known_compounds:
                                reaction_order(reaction_product, reaction_ordered, reactions_compounds, alreayd_known_compounds)

            for input_compound in pathway_input_compounds:
                reaction_ordered = []
                alreayd_known_compounds = []
                reaction_order(input_compound, reaction_ordered, reactions_compounds, alreayd_known_compounds)
                if reaction_ordered != []:
                    pathway_node.misc["REACTIONS-ORDER"] = [','.join(reaction_ordered)]

        try:
            types = dict_values["TYPES"]
            _setType(types, pathway_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, pathway_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, pathway_id, padmet)
        except KeyError:
            pass
        try:
            subPathways = dict_values["IN-PATHWAY"]
            for subPathway in subPathways:
                #add the hierachization info, current pathway is_in_pathway subpathway
                is_in_pathway_rlt = Relation(pathway_id, "is_in_pathway", subPathway)
                list_of_relation.append(is_in_pathway_rlt)
        except KeyError:
            pass
        try:
            subNodes = dict_values["REACTION-LIST"]
            for subNode in subNodes:
                #add the hierachization info, Reaction/pathway is_in_pathway current pathway
                is_in_pathway_rlt = Relation(subNode, "is_in_pathway", pathway_id)
                list_of_relation.append(is_in_pathway_rlt)
        except KeyError:
            pass
    if verbose: print("")

def compounds_parser(filePath, padmet, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to compounds.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """        
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","INCHI-KEY","MOLECULAR-WEIGHT","SMILES",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

    compounds = []
    count = 0
    nb_cpds = str(len(list(dict_data.keys())))
    for compound_id, dict_values in dict_data.items():
        count += 1
        if verbose:
            print("\r%s/%s: %s" %(count, nb_cpds, compound_id),end="", flush=True)
        compound_node = Node("compound", compound_id)
        padmet.dicOfNode[compound_id] = compound_node
        compounds.append(compound_id)
        try:
            compound_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            compound_node.misc["INCHI-KEY"] = dict_values["INCHI-KEY"]
        except KeyError:
            pass
        try:
            compound_node.misc["MOLECULAR-WEIGHT"] = dict_values["MOLECULAR-WEIGHT"]
        except KeyError:
            pass
        try:
            compound_node.misc["SMILES"] = dict_values["SMILES"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types,compound_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, compound_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, compound_id, padmet)
        except KeyError:
            pass
    if verbose: print("")

    return compounds

def rnas_parser(filePath, padmet, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to compounds.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME", "TYPES", "GENE", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

    count = 0
    rna_genes = {}
    nb_rnas = str(len(list(dict_data.keys())))
    for rna_id, dict_values in dict_data.items():
        count += 1
        if verbose:
            print("\r%s/%s: %s" %(count, nb_rnas, rna_id),end="", flush=True)
        try:
            rna_genes[rna_id] = dict_values["GENE"]
        except KeyError:
            rna_genes[rna_id] = None

    if verbose: print("")

    return rna_genes

def genes_parser(filePath, padmet, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to genes.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """        
    dict_data = {}
    #k='ACCESSION-1', v ='PRODUCT'
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","ACCESSION-1","CENTISOME-POSITION","LEFT-END-POSITION","RIGHT-END-POSITION","SYNONYMS","TRANSCRIPTION-DIRECTION","PRODUCT"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

    count = 0
    nb_genes = str(len(list(dict_data.keys())))
    map_gene_ids = {}
    for current_id, dict_values in dict_data.items():
        try:
            gene_id = dict_values.get("ACCESSION-1",[current_id])[0]
            map_gene_ids[current_id] = gene_id
            count += 1
            if verbose: 
                print("\r%s/%s: %s" %(count, nb_genes, gene_id), end="", flush=True)

            gene_node = Node("gene", gene_id)
            padmet.dicOfNode[gene_id] = gene_node
            try:
                if dict_values["COMMON-NAME"][0] != gene_id:
                    gene_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
            except KeyError:
                pass
            try:
                if dict_values["TRANSCRIPTION-DIRECTION"][0] == "-":
                    gene_node.misc["TRANSCRIPTION-DIRECTION"] = ["NEGATIVE"]
                elif dict_values["TRANSCRIPTION-DIRECTION"][0] == "+":
                    gene_node.misc["TRANSCRIPTION-DIRECTION"] = ["POSITIVE"]
            except KeyError:
                pass
            try:
                gene_node.misc["CENTISOME-POSITION"] = dict_values["CENTISOME-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["LEFT-END-POSITION"] = dict_values["LEFT-END-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["RIGHT-END-POSITION"] = dict_values["RIGHT-END-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["RIGHT-END-POSITION"] = dict_values["RIGHT-END-POSITION"]
            except KeyError:
                pass
            try:
                syns = dict_values["SYNONYMS"]
                _setSyns(syns, gene_id, padmet)
            except KeyError:
                pass
        except KeyError:
            pass
    if verbose: print("")

    return map_gene_ids

def proteins_parser(filePath, padmet, compounds, rnas, id_classes, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to proteins.dat
    padmet: padmet.PadmetRef
        padmet instance
    compounds: list
        list of known compounds in the species
    verbose: bool
        if True print information
    """        
    dict_data = {}
    dict_protein_component_id = {}
    dict_protein_gene_id = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","INCHI-KEY","MOLECULAR-WEIGHT","SMILES",\
                "TYPES", "SPECIES", "SYNONYMS", "DBLINKS", "GENE", "GO-TERMS", "COMPONENTS","COMPONENT-OF"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            
    count = 0
    nb_proteins = str(len(list(dict_data.keys())))
    for protein_id, dict_values in dict_data.items():
        count += 1
        if verbose:
            print("\r%s/%s" %(count, nb_proteins), end="", flush=True)

        try:
            dict_protein_component_id[protein_id] = dict_values["COMPONENTS"]
        except KeyError:
            pass
        try:
            dict_protein_gene_id[protein_id] = dict_values["GENE"]
        except KeyError:
            pass

    def get_gene_id(protein_id, list_of_components, compounds, rnas, id_classes):
        """
        Get the gene ID corresponding to the protein or the complex.
        Some complex contain other complex, so we use recursivity to get the complex required.
        Also some compounds or RNAs can be part of a complex so we take only component in the gene ID lists.
        """
        genes_associated = set()
        for component in list_of_components:
            if component in rnas:
                if rnas[component]:
                    genes_associated.update(rnas[component])

            elif component not in compounds and component not in id_classes:
                try:
                    genes_associated.update(dict_protein_gene_id[component])
                except KeyError:
                    if component in dict_protein_component_id:
                        get_gene_id(component, dict_protein_component_id[component], compounds, rnas, id_classes)
                        genes_associated.update(dict_protein_gene_id[component])
        dict_protein_gene_id[protein_id] = genes_associated

    # Extract protein and complexes.
    for protein_id, list_of_components in dict_protein_component_id.items():
        get_gene_id(protein_id, list_of_components, compounds, rnas, id_classes)

    if verbose: print("")

    return dict_protein_gene_id

def enzrxns_parser(filePath, padmet, dict_protein_gene_id, source, verbose = False):
    """
    Parameters
    ----------
    filePath: str
        path to enzrxns.dat
    padmet: padmet.PadmetRef
        padmet instance
    verbose: bool
        if True print information
    """        
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            if len(line.split(" - ")) > 1:
                attrib = line.split(" - ")[0]
                value = " - ".join(line.split(" - ")[1:])
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","ENZYME","REACTION","BASIS-FOR-ASSIGNMENT"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

    count = 0
    nb_enzrxns = str(len(list(dict_data.keys())))
    for current_id, dict_values in dict_data.items():
        count += 1

        if verbose:
            print("\r%s/%s" %(count, nb_enzrxns), end="", flush=True)

        rxn_id = dict_values["REACTION"][0]
        names = dict_values.get("COMMON-NAME",[])
        for name in names: 
            if name.endswith("_"):
                names[names.index(name)] = name[:-1]
        try:
            protein = dict_values["ENZYME"][0]
        except KeyError:
            pass
        try:
            rxn_node = padmet.dicOfNode[rxn_id]
            try:
                [rxn_node.misc["COMMON-NAME"].append(name) for name in names if name not in rxn_node.misc["COMMON-NAME"]]
            except KeyError:
                rxn_node.misc["COMMON-NAME"] = names
            try:
                genes_id = dict_protein_gene_id[protein]
                try:
                    assignment = dict_values["BASIS-FOR-ASSIGNMENT"][0]
                    if assignment.startswith(":"): assignment = assignment[1:]
                except KeyError:
                    assignment = "NA"
                for gene_id in genes_id:
                    is_linked_rlt = Relation(rxn_id, "is_linked_to", gene_id, {"SOURCE:ASSIGNMENT":[source+":"+assignment]})
                    list_of_relation.append(is_linked_rlt)
            except KeyError:
                pass
        except KeyError:
            pass
    if verbose: print("")


def _setType(types, current_id, padmet):
    """
    For id current_id, create relation 'is_a_class' to each id in types
    add new nodes in padmet and relation in list_of_relation
    Parameters
    ----------
    types: list
        list of classe id associated to id 'current_id'
    current_id: str
        element id to link to all classe id in types
    padmet: padmet.PadmetRef
        padmet instance
    """        
    for subClass_id in types:
        #Type allow the hierachization of the current
        #XX is_a_class type
        #if type is not already in padmet, create a new class node (subClass_node)
        #subClass_node id == type
        #create a relation xx is_a_class type
        try:
            subClass_node = padmet.dicOfNode[subClass_id]
        except KeyError:
            subClass_node = Node("class", subClass_id)
            padmet.dicOfNode[subClass_id] = subClass_node
        is_a_class_rlt = Relation(current_id, "is_a_class", subClass_id)
        list_of_relation.append(is_a_class_rlt)
    
def _setSyns(syns, current_id, padmet):
    """
    For id current_id, create relation 'has_name' to a node name 'current_id'_names
    store a list of synonymous from syns list in the node name 'current_id'_names
    add new node in padmet and relation in list_of_relation

    Parameters
    ----------
    syns: list
        list of synonymous of current_id
    current_id: str
        element id to link to node name 'current_id'_names
    padmet: padmet.PadmetRef
        padmet instance
    """        
    name_id = current_id+"_names"
    try:
        name_node = padmet.dicOfNode[name_id]
    except KeyError:
        #create node name
        name_node = Node("name", name_id, {"LABEL":[]})
        padmet.dicOfNode[name_id] = name_node
        has_name_rlt = Relation(current_id, "has_name", name_id)
        list_of_relation.append(has_name_rlt)
    [name_node.misc["LABEL"].append(syn) for syn in syns 
    if syn not in name_node.misc["LABEL"]]

def _setXrefs(xrefs, current_id, padmet):
    """
    For id current_id, create relation 'has_xref' to a node xref 'current_id'_xrefs
    store a list of external reference from xrefs list in the node xref 'current_id'_xrefs
    parse each external reference with regex_xref
    add new node in padmet and relation in list_of_relation

    Parameters
    ----------
    xrefs: list
        list of external reference of current_id
    current_id: str
        element id to link to node xref 'current_id'_xref
    padmet: padmet.PadmetRef
        padmet instance
    """        
    xref_id = current_id+"_xrefs"
    try:
        xref_node = padmet.dicOfNode[xref_id]
    except KeyError:
        #create node xref
        xref_node = Node("xref", xref_id)
        padmet.dicOfNode[xref_id] = xref_node
        has_xref_rlt = Relation(current_id, "has_xref", xref_id)
        list_of_relation.append(has_xref_rlt)

    for xref in xrefs:
        #an xref is like: (REFSEQ "NP_417401" NIL NIL NIL NIL NIL)
        #in this example DB = REFSEQ and ID = NP_417401
        #update node xref, with in misc k = DB and v = [ID]
        #node id is created by incrementing meta_max_id
        xref_search = regex_xref.search(xref)
        if xref_search is not None:
            xref_dict = xref_search.groupdict()
            db = xref_dict["DB"]
            _id = xref_dict["ID"]
        else:
            db = "GO-TERMS"
            _id = xref
        if db in list(xref_node.misc.keys()) and _id not in xref_node.misc[db]:
            xref_node.misc[db].append(_id)
        else:
            xref_node.misc[db] = [_id]

def enhance_db(metabolic_reactions, padmet, with_genes, verbose = False):
    """
    Parse sbml metabolic_reactions and add reactions in padmet
    if with_genes: add also genes information

    Parameters
    ----------
    metabolic_reactions: str
        path to sbml metabolic-reactions.xml
    padmet: padmet.PadmetRef
        padmet instance
    with_genes: bool
        if true alos add genes information.

    Returns
    -------
    padmet.padmetRef:
        padmet instance with pgdb within pgdb + metabolic-reactions.xml data
    """        
    
    print("loading sbml file: %s" %metabolic_reactions)
    reader = libsbml.SBMLReader()
    document = reader.readSBML(metabolic_reactions)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    #recovere the reactions that are not in the basic metacyc but in the sbml file
    #use the reactions_name instead of ids because the ids are encoded, the name is the non-encoded version of the id
    padmet_reactions_id = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"])
    reaction_to_add = [reaction for reaction in listOfReactions 
    if reaction.getName() not in padmet_reactions_id]
    count = 0
    if verbose: print(str(len(reaction_to_add))+" reactions to add")
    for reactionSBML in reaction_to_add:
        count += 1
        reaction_id = reactionSBML.getName()
        if verbose: print(str(count)+"/"+str(len(reaction_to_add))+"\t"+reaction_id)
        if reactionSBML.getReversible():
            reaction_dir = "REVERSIBLE"
        else:
            reaction_dir = "LEFT-TO-RIGHT"
        try:
            reaction_node = padmet.dicOfNode[reaction_id]
        except KeyError:
            reaction_node = Node("reaction", reaction_id, {"DIRECTION": [reaction_dir]})
            padmet.dicOfNode[reaction_id] = reaction_node
        reactants = reactionSBML.getListOfReactants()
        for reactant in reactants: #convert ids
            reactant_id, _type, reactant_compart = sbmlPlugin.convert_from_coded_id(reactant.getSpecies())
            if reactant_id not in list(padmet.dicOfNode.keys()):
                reactant_node = Node("compound",reactant_id)
                padmet.dicOfNode[reaction_id] = reactant_node
            reactant_stoich = reactant.getStoichiometry()
            consumes_rlt = Relation(reaction_id,"consumes",reactant_id, {"STOICHIOMETRY":[reactant_stoich], "COMPARTMENT": [reactant_compart]})
            list_of_relation.append(consumes_rlt)

        products = reactionSBML.getListOfProducts()
        for product in products:
            product_id, _type, product_compart = sbmlPlugin.convert_from_coded_id(product.getSpecies())
            if product_id not in list(padmet.dicOfNode.keys()):
                product_node = Node("compound",product_id)
                padmet.dicOfNode[product_id] = product_node
            product_stoich = product.getStoichiometry()
            produces_rlt = Relation(reaction_id,"produces",product_id,{"STOICHIOMETRY": [product_stoich], "COMPARTMENT": [product_compart]})
            list_of_relation.append(produces_rlt)
        
        if with_genes:
            notes = sbmlPlugin.parseNotes(reactionSBML)
            if "GENE_ASSOCIATION" in list(notes.keys()):
                #Using sbmlPlugin to recover all genes associated to the reaction
                listOfGenes = sbmlPlugin.parseGeneAssoc(notes["GENE_ASSOCIATION"][0])
                if len(listOfGenes) != 0:
                    for gene in listOfGenes:
                        try:
                            #check if gene already in the padmet
                            padmet.dicOfNode[gene]
                        except TypeError:
                            gene_node = Node("gene",gene)
                            padmet.dicOfNode[gene] = gene_node
                        is_linked_rlt = Relation(reaction_id, "is_linked_to", gene)
                        list_of_relation.append(is_linked_rlt)
    return padmet

def map_gene_id(dict_protein_gene_id, map_gene_ids):
    """
    Map gene ID created by Pathway Tools with gene ID from the data.
    Automatically Pathway Tools uppercased all the letter in gene ID.
    So we need to do this mapping to retrieve the unuppercased gene ID.
    """
    mapped_dict_protein_gene_id = {}

    for prot_id in dict_protein_gene_id:
        mapped_gene_ids = set()
        for gene_id in dict_protein_gene_id[prot_id]:
            if gene_id in map_gene_ids:
                mapped_gene_ids.add(map_gene_ids[gene_id])
        mapped_dict_protein_gene_id[prot_id] = mapped_gene_ids

    return mapped_dict_protein_gene_id

