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
from padmet.classes import PadmetSpec, PadmetRef,instantiate_padmet
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
        padmet_to_update = instantiate_padmet("PadmetRef", None, db, version, verbose)

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
        padmet_to_update = instantiate_padmet("PadmetSpec", padmetRef_file, db, version, verbose)

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
