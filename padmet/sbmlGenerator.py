#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of padmet.

padmet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
The module sbmlGenerator contains functions to generate sbml files from padmet and txt
usign the libsbml package
"""
from padmetSpec import PadmetSpec
from padmetRef import PadmetRef
import sbmlPlugin
try:
    from libsbml import *
except:
    print("package libsbml needed, use this cmd:\n pip install "
          + "python-libsbml")
    exit()

def check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value == None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + OperationReturnValue_toString(value).strip() + '"'
            raise SystemExit(err_msg)
    else:
        return

  
def padmet_to_sbml(padmet_file, output, obj_fct = None, sbml_lvl = 2, sbml_version = 1, verbose = False):
    """
    Convert padmet file to sbml file.
    Specificity: 
    - ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param padmet_file: the pathname to the padmet file to convert
    @param output: the pathname to the sbml file to create
    @param obj_fct: the identifier of the objection function, the reaction to test in FBA
    @param sbml_lvl: the sbml level
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type padmet_file, output, verbose: str
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    #create an empty sbml model
    document = SBMLDocument(sbml_lvl, sbml_version)
    model = document.createModel()
    #load padmet
    padmet = PadmetSpec(padmet_file)

    BOUNDARY_ID = 'C-BOUNDARY'
    default_lower_bound = -1000
    default_upper_bound = 1000
    math_ast = parseL3Formula('FLUX_VALUE')
    check(math_ast, 'create AST for rate expression')

    # Create a unit definition
    mmol_per_gDW_per_hr = model.createUnitDefinition()
    check(mmol_per_gDW_per_hr, 'create unit definition')
    check(mmol_per_gDW_per_hr.setId('mmol_per_gDW_per_hr'), 'set unit definition id')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create mole unit')
    check(unit.setKind(UNIT_KIND_MOLE), 'set unit kind')
    check(unit.setScale(-3), 'set unit scale')
    check(unit.setMultiplier(1), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create gram unit')
    check(unit.setKind(UNIT_KIND_GRAM), 'set unit kind')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setMultiplier(1), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create second unit')
    check(unit.setKind(UNIT_KIND_SECOND), 'set unit kind')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setMultiplier(0.00027777), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')


    #generator of tuple: (x,y) x=species id,y=value of compart, if not defined=""
    species = [(rlt.id_out, rlt.misc.get("COMPARTMENT",[None])[0]) for rlt in padmet.getAllRelation() 
    if rlt.type in ["consumes","produces"]]
    if verbose: print("%s species" %len(species))
    #compart_dict: k = id_encoded, v = original id
    compart_dict = {}    
    #species_dict: k = species_id_encoded, v = dict: k' = {species_id, compart, name}, v' = value or None 
    species_dict = {}
    for species_id, compart in species:
        #encode id for sbml
        species_id_encoded = sbmlPlugin.convert_to_coded_id(species_id, "M", compart)
        
        #encode compart id for sbml
        #try to get the common_name, if non value return None
        name = padmet.dicOfNode[species_id].misc.get("COMMON_NAME",[species_id])[0]
        #update dicts
        species_dict[species_id_encoded] = {"species_id":species_id, "compart":compart, "name":name}
        
    
    for k, v in species_dict.iteritems():
        compart = v["compart"]
        name = v["name"]
        
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(k), 'set species id %s' %k)
        check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
        #check(s.setMetaId(metaId), 'set species MetaId %s' %metaId)
        if name is not None:
            check(s.setName(name), 'set species Name %s' %name)
        if compart is not None:
            compart_encoded = sbmlPlugin.convert_to_coded_id(compart)
            compart_dict[compart_encoded] = compart
            check(s.setCompartment(compart_encoded), 'set species compartment %s' %compart_encoded)
            if compart == BOUNDARY_ID:
                check(s.setBoundaryCondition(True), 'set boundaryCondition to True')

    for k, v in compart_dict.iteritems():
        compart = model.createCompartment()
        check(compart,'create compartment')
        check(compart.setId(k),'set compartment id %s' %k)
        if v == "c":
            check(compart.setName("cytosol"),'set compartment name cytosol')
        elif v == "e":
            check(compart.setName("extracellular"),'set compartment name extracellular')
        elif v == "p":
            check(compart.setName("periplasm"),'set compartment name periplasm')
        elif v != k:
            check(compart.setName(v),'set compartment id %s' %v)


    if obj_fct is not None:
        obj_fct_encoded = sbmlPlugin.convert_to_coded_id(obj_fct)
        if verbose: print("the objectif reaction is: %s initialy: %s" %(obj_fct_encoded, obj_fct))
        
    reactions = [node for node in padmet.dicOfNode.itervalues() if node.type == "reaction"]
    nb_reactions = str(len(reactions))    
    # Create reactions
    if verbose: print("%s reactions" %nb_reactions)

    for rNode in reactions:
        rId = rNode.id
        rId_encoded = sbmlPlugin.convert_to_coded_id(rId,"R")
        rName = rNode.misc.get("COMMON_NAME",[rId])[0]

        #generator of tuple (reactant_id,stoichiometry,compart)
        try:
            consumed = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
            for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "consumes")
        except TypeError:
            print rId
            exit()
            
        #generator of tuple (product_id,stoichiometry,compart)        
        produced = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
        for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "produces")
        direction = rNode.misc["DIRECTION"][0]
        if direction == "LEFT-TO-RIGHT":
            reversible = False
        #include if direction = unknown
        else:
            reversible = True
        
        reaction = model.createReaction()
        check(reaction, 'create reaction')
        check(reaction.setId(rId_encoded), 'set reaction id %s' %rId_encoded)
        if rName is not None:
            check(reaction.setName(rName), 'set reaction name %s' %rName)
        check(reaction.setReversible(reversible), 'set reaction reversibility flag %s' %reversible)
        
        kinetic_law = reaction.createKineticLaw()
        check(kinetic_law, 'create kinetic law')
        check(kinetic_law.setMath(math_ast), 'set math on kinetic law')
        #add parameter flux_value
        flux_value_k = kinetic_law.createParameter()
        check(flux_value_k, 'create parameter flux_value_k')
        check(flux_value_k.setId('FLUX_VALUE'), 'set parameter flux_value_k id')
        check(flux_value_k.setValue(0), 'set parameter flux_value_k value')
        check(flux_value_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter flux_value_k units')
        #add parameter upper/lower_bound, lower value depend on reversibility
        upper_bound_k = kinetic_law.createParameter()
        check(upper_bound_k, 'create parameter upper_bound_k')
        check(upper_bound_k.setId('UPPER_BOUND'), 'set parameter upper_bound_k')
        check(upper_bound_k.setValue(default_upper_bound),'set parameter upper_bounp_k value')
        check(upper_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter uppper_bound_k units')

        if reversible:
            lower_bound_k = kinetic_law.createParameter()
            check(lower_bound_k, 'create parameter lower_bound_k')
            check(lower_bound_k.setId('LOWER_BOUND'), 'set parameter lower_bound_k id')
            check(lower_bound_k.setValue(default_lower_bound), 'set parameter lower_bound_k value')
            check(lower_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter lower_bound_k units')
        else:
            lower_bound_k = kinetic_law.createParameter()
            check(lower_bound_k, 'create parameter lower_bound_k')
            check(lower_bound_k.setId('LOWER_BOUND'), 'set parameter lower_bound_k id')
            check(lower_bound_k.setValue(0), 'set parameter lower_bound_k value')
            check(lower_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter lower_bound_k units')
        #objective_coeeficient
        if rId == obj_fct:
            obj_fct_k = kinetic_law.createParameter()
            check(obj_fct_k, 'create parameter obj_fct_k')
            check(obj_fct_k.setId('OBJECTIVE_COEFFICIENT'), 'set parameter obj_fct_k id')
            check(obj_fct_k.setValue(1), 'set parameter obj_fct_k value')
        else:
            obj_fct_k = kinetic_law.createParameter()
            check(obj_fct_k, 'create parameter obj_fct_k')
            check(obj_fct_k.setId('OBJECTIVE_COEFFICIENT'), 'set parameter obj_fct_k id')
            check(obj_fct_k.setValue(0), 'set parameter obj_fct_k value')

        for cId, stoich, compart in consumed:
            cId_encoded = sbmlPlugin.convert_to_coded_id(cId,"M",compart)
            try:
                stoich = float(stoich)
            #for case stoich = n
            except ValueError:
                stoich = float(1)
            species_ref = reaction.createReactant()
            check(species_ref, 'create reactant')
            check(species_ref.setSpecies(cId_encoded), 'assign reactant species %s' %cId_encoded)
            check(species_ref.setStoichiometry(stoich), 'set stoichiometry %s' %stoich)

        for pId, stoich, compart in produced:
            pId_encoded = sbmlPlugin.convert_to_coded_id(pId,"M",compart)
            try:
                stoich = float(stoich)
            except ValueError:
                stoich = float(1)
            species_ref = reaction.createProduct()
            check(species_ref, 'create product')
            check(species_ref.setSpecies(pId_encoded), 'assign product species %s' %pId_encoded)
            check(species_ref.setStoichiometry(stoich), 'set stoichiometry %s' %stoich)

        try:
            linked_genes = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rId, [])
            if rlt.type == "is_linked_to"])
            if len(linked_genes) != 0:
                if verbose: print rId, "is linked to", len(linked_genes), "genes."
                #Try to recover linked genes expressions from suppData
                linked_genes_from_suppData = " or ".join([padmet.dicOfNode[rlt.id_out].misc.get("GENE_ASSOCIATION", [])[0] 
                for rlt in padmet.dicOfRelationIn.get(rId, [])
                if rlt.type == "has_suppData"])
                
                notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
                if linked_genes_from_suppData:
                    notes += "<p>"+"GENE_ASSOCIATION:" + linked_genes_from_suppData + "</p>"
                else:
                    linked_genes = " or ".join(linked_genes)
                    notes += "<p>"+"GENE_ASSOCIATION:" + linked_genes + "</p>"

                notes += "</body>"
                check(reaction.setNotes(notes), 'set notes %s' %notes)
        except IndexError:
            pass

    if verbose: print("Done, creating sbml file: %s" %output)
    return writeSBMLToFile(document, output)   

#################################

def reactions_to_SBML(reactions_file, output, padmetRef_file, verbose = False):
    """
    convert a list of reactions to sbml format based on a given padmet of reference.
    - ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param reactions_file:the pathname to the file containing the reactions ids, 1/line
    @param padmetRef_file: the pathname to the file padmet of reference
    @param output: the pathname to the sbml file to create
    @param sbml_lvl: the sbml level
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type reactions_file, output, padmetRef_file, verbose: str
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    padmetRef = PadmetRef(padmetRef_file)
    with open(reactions_file,'r') as f:
        set_of_reactions_id = set(f.read().splitlines())
    
    all_rxn = set([k for k,v in padmetRef.dicOfNode.iteritems() if v.type == "reaction"])
    rxn_not_in_ref = set_of_reactions_id.difference(all_rxn)
    if len(rxn_not_in_ref) == len(set_of_reactions_id):
        raise SystemError("None of the reactions are in padmetRef")
    else:
        for rxn_id in rxn_not_in_ref:
            if verbose: print("%s not in padmetRef" % rxn_id)
            set_of_reactions_id.pop(rxn_id)

    diff_rxn = all_rxn.difference(set_of_reactions_id)
    for rxn_id in diff_rxn:
        padmetRef.dicOfNode.pop(rxn_id)
    padmet_to_SBML(padmetRef, output, verbose = verbose)
    

    return writeSBMLToFile(document, output)

def compounds_to_sbml(compounds_file, output, padmetRef_file = None, padmetSpec_file = None, sbml_lvl = 2, sbml_version = 1, verbose = False):
    """
    convert a list of compounds to sbml format
    if compart_name is not None, then the compounds id will by: M_originalID_compart_name
    if verbose and specified padmetRef and/or padmetSpec: will check if compounds are in one of the padmet files
    Ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param compounds_file: the pathname to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
    @param output: the pathname to the sbml file to create
    @param padmetRef_file: the pathname to the file padmet of reference
    @param padmetRef_file: the pathname to the file padmet of a species
    @param compart_name: the default compart to concatenate
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type compounds_file, output, padmetRef_file, padmetSpec_file, verbose: str
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    if verbose:
        if padmetRef_file is not None:
            padmetRef = PadmetRef(padmetRef_file)
        if padmetSpec_file is not None:
            padmetSpec = PadmetSpec(padmetSpec_file)
    document = SBMLDocument(sbml_lvl, sbml_version)
    model = document.createModel()

    with open(compounds_file,'r') as f:
        dict_of_compounds = dict([line.split("\t") for line in f.read().splitlines()])
    all_compart = [sbmlPlugin.convert_to_coded_id(compart) for compart in set(dict_of_compounds.values())]
    for k in all_compart:
        compart = model.createCompartment()
        check(compart,'create compartment')
        check(compart.setId(k),'set compartment id %s' %k)

    if verbose:
        in_padmetRef, in_padmetSpec = set(), set()
        if padmetRef_file is not None:
            for sId in dict_of_compounds.iterkeys():
                try:
                    padmetRef.dicOfNode[sId]
                    in_padmetRef.add(sId)
                except KeyError:
                    pass
        if padmetSpec_file is not None:
            for sId in dict_of_compounds.iterkeys():
                try:
                    padmetSpec.dicOfNode[sId]
                    in_padmetSpec.add(sId)
                except KeyError:
                    pass
        in_both = in_padmetRef.intersection(in_padmetSpec)
        if verbose:
            for sId in dict_of_compounds.iterkeys():
                if sId in in_both and (padmetRef_file is not None and padmetSpec_file is not None):
                    print(sId+" in padmetRef and padmetSpec")
                elif sId in in_padmetRef and padmetRef_file is not None:
                    print(sId+" in padmetRef only")
                elif sId in in_padmetSpec and padmetSpec_file is not None:
                    print(sId+" in padmetSpec only")
                else:
                    if padmetRef_file is not None or padmetSpec_file is not None:
                        print(sId+" not in padmetRef or padmetSpec")
              
    nb_species = str(len(dict_of_compounds.keys()))
    count = 0
    for sName, sCompart in dict_of_compounds.iteritems():
        sId_encoded = sbmlPlugin.convert_to_coded_id(sName,"M",sCompart)
        count += 1
        if verbose: print(str(count)+"/"+nb_species+"\t"+sName)
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(sId_encoded), 'set species id')
        check(s.setName(sName), 'set species name')
        check(s.setCompartment(sbmlPlugin.convert_to_coded_id(sCompart)), 'set species compartment')
    
    return writeSBMLToFile(document, output)
