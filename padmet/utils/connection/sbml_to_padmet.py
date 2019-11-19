# -*- coding: utf-8 -*-
"""
Description:
    There are 3 cases of convertion sbml to padmet:
    
    1./ Creation of a reference database in padmet format from sbml(s) (or updating one with new(s) sbml(s))
    First usage, padmetRef is the padmetRef to create or to update. If it's an update case, the output
    can be used to create a new padmet, if output None, will overwritte the input padmetRef.
    
    2./ Creation of a padmet representing an organism in padmet format from sbml(s) (or updating one with new(s) sbml(s))
    2.A/ Without a database of reference:
    Second usage, padmetSpec is the padmetSpec to create or update. If it's an update case, the output
    can be used to create a new padmet, if output None, will overwritte the input padmetSpec.
    
    2.B/ With a database of refence:
    Third usage, padmetSpec is the padmetSpec to create or update. If it's an update case, the output
    can be used to create a new padmet, if output None, will overwritte the input padmetSpec.
    padmetRef is the padmet representing the database of reference.
    
    It is possible to define a specific policy and info for the padmet. To learn more about
    policy and info check doc of lib.padmetRef/Spec.
    if the ids of reactions/compounds are not the same between padmetRef and the sbml, it is possible to use
    a dictionnary of association (sbml_id padmetRef_id)
    with one line = 'id_sbml \t id_padmetRef'
    Finally if a reaction from sbml is not in padmetRef, it is possible to force the copy and creating
    a new reaction in padmetSpec with the arg -f

"""
from padmet.classes import PadmetSpec, PadmetRef
from datetime import datetime
import os
        
def sbml_to_padmetRef(sbml, padmetRef_file, output=None, db="NA", version="NA", verbose=False):
    """
    if padmetRef, not padmetSpec:
        if padmetRef exist, instance PadmetRef
        else init PadmetRef
        update padmetRef
    if padmetSpec:
        if padmetRef, check if exist else raise Error
        if padmetSpec exist, instance PadmetSpec
        else init PadmetSpec
        update padmetSpec using padmetRef if padmetRef
    
    #TODO
    """
    if output is None:
        output = padmetRef_file
    if os.path.isdir(sbml):
        sbml_files = [os.path.join(sbml,_f) for _f in next(os.walk(sbml))[2] if _f.endswith(".sbml") or _f.endswith(".xml")]
    else:
        sbml_files = sbml.split(";")

    if os.path.isfile(padmetRef_file):
        padmet_to_update = PadmetRef(padmetRef_file)
    else:
        padmet_to_update = create_padmet_instance(padmetRef_file, "PadmetRef", db, version)

    for sbml_file in sbml_files:
        if verbose:
            print("Updating padmet from %s" %os.path.basename(sbml_file))
        padmet_to_update.updateFromSbml(sbml_file, verbose)

    padmet_to_update.generateFile(output)

def sbml_to_padmetSpec(sbml, padmetSpec_file, padmetRef_file=None, output=None, mapping=None, mapping_tag="_dict.csv", source_tool=None, source_category=None, db="NA", version="NA", verbose=False):
    """
    Convert 1 - n sbml to padmet file.
    sbml var is file or dir
    padmetSpec_file: path to new padmet file to create or old padmet to update
    padmetRef_file: path to database of reference to use for data standardization
    output: path to new padmet file, if none, overwritte padmetSpec_file
    source_tool: tool used to create this sbml(s) ex Orthofinder
    source_category: Category of the tool ex: orthology
    if new padmet without padmetRef:
        db: database used ex: metacyc, bigg
        version: version of the database, 23, 18...
    


    if padmetRef, not padmetSpec:
        if padmetRef exist, instance PadmetRef
        else init PadmetRef
        update padmetRef
    if padmetSpec:
        if padmetRef, check if exist else raise Error
        if padmetSpec exist, instance PadmetSpec
        else init PadmetSpec
        update padmetSpec using padmetRef if padmetRef
    
    #TODO
    """
    if output is None:
        output = padmetSpec_file
    #if sbml is a dir: sbml_files are all files with extension .sbml or .xml within dir
    #else: sbml = my_sbml.sbml or sbml= my_sbml1.sbml;my_sbml2.sml
    if os.path.isdir(sbml):
        sbml_files = [os.path.join(sbml,_f) for _f in next(os.walk(sbml))[2] if _f.endswith(".sbml") or _f.endswith(".xml")]
    else:
        sbml_files = [sbml]
    
    #PadmetRef used for mapping and data standardization
    if padmetRef_file:
        padmetRef = PadmetRef(padmetRef_file)
    else:
        padmetRef = None

    if os.path.isfile(padmetSpec_file):
        padmet_to_update = PadmetSpec(padmetSpec_file)
    else:
        padmet_to_update = create_padmet_instance(padmetSpec_file, "PadmetSpec", db, version, padmetRef)

    #if sbml is a directory, recover all file path in a list. if no => only one file: create a list with only this file
    #sbml_mapping_dict = {'/path/to/my_sbml1.sbml': '/path/to/my_sbml1_dict.csv' // None}  
    sbml_mapping_dict = {}
    if os.path.isdir(sbml):
        for sbml_file in sbml_files:
            mapping_file = os.path.splitext(sbml_file)[0] + mapping_tag
            if not os.path.isfile(mapping_file) or not padmetRef:
                mapping_file = None
            sbml_mapping_dict[sbml_file] = mapping_file
    else:
        sbml_mapping_dict[sbml] = mapping

    if len(list(sbml_mapping_dict.keys())) == 0:
        raise IOError("No sbml found based on %s" %sbml)

    for sbml_file, mapping_file in list(sbml_mapping_dict.items()):
        if mapping_file:
            force = False
        else:
            force = True

        if verbose:
            if mapping_file:
                print("Updating %s from %s using mapping dictionnary %s" %(os.path.basename(padmetSpec_file),os.path.basename(sbml_file),os.path.basename(mapping_file)))
            else:
                print("Updating %s from %s" %(os.path.basename(padmetSpec_file),os.path.basename(sbml_file)))

        padmet_to_update.updateFromSbml(sbml_file=sbml_file, padmetRef=padmetRef, mapping_file=mapping_file, verbose=verbose, force=force, source_category=source_category, source_tool=source_tool)

    padmet_to_update.generateFile(output)



def create_padmet_instance(padmet_file, padmet_type, db, version, padmetRef=None):
    """
    #TODO
    """
    if padmet_type not in ["PadmetRef","PadmetSpec"]:
        raise TypeError('padmet_type must be in ["PadmetRef","PadmetSpec"], given:%s' %padmet_type)
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")

    if padmet_type == "PadmetSpec":
        padmet = PadmetSpec()
    elif padmet_type == "PadmetRef":
        padmet = PadmetRef()
        
    if padmetRef:
        padmet.setInfo(padmetRef)
        padmet.info["PADMET"]["creation"] = today_date
        padmet.setPolicy(padmetRef)
    else:
        POLICY_IN_ARRAY = [['class','is_a_class','class'], ['class','has_name','name'], ['class','has_xref','xref'], ['class','has_suppData','suppData'],
                        ['compound','is_a_class','class'], ['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                        ['gene','is_a_class','class'], ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                        ['pathway','is_a_class','class'], ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'],
                        ['protein','is_a_class','class'], ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                        ['protein','is_in_species','class'],
                        ['reaction','is_a_class','class'], ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],
                        ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'],
                        ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'],
                        ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'],
                        ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
        dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":db,"version":version}}
        padmet.setInfo(dbNotes)
        padmet.setPolicy(POLICY_IN_ARRAY)
    return padmet

