# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in two sbml.

    Returns if a reaction is missing

    And if a reaction with the same id is using different species or different reversibility

"""
from cobra.io.sbml import create_cobra_model_from_sbml_file

def compare_sbml(sbml1_path, sbml2_path):
    """
    Compare 2 sbml, print nb of metabolites and reactions.
    If reaction missing print reaction id, and reaction formula.

    Parameters
    ----------
    sbml1_path: str
        path to the first sbml file to compare
    sbml2_path: str
        path to the second sbml file to compare
    """
    sbml_1 = create_cobra_model_from_sbml_file(sbml1_path)
    sbml_2 = create_cobra_model_from_sbml_file(sbml2_path)
    
    print("sbml1:")
    print("metabolites: %s" %(len(sbml_1.metabolites)))
    print("reactions: %s" %(len(sbml_1.reactions)))
    print("sbml2:")
    print("metabolites: %s" %(len(sbml_2.metabolites)))
    print("reactions: %s" %(len(sbml_2.reactions)))
    not_in1 = [i for i in sbml_2.reactions if i not in sbml_1.reactions]
    print("reactions not in sbml1: %s" %len(not_in1))
    for i in not_in1:
        print("\t%s" %i)
    not_in2 = [i for i in sbml_1.reactions if i not in sbml_2.reactions]
    print("reactions not in sbml2: %s" %len(not_in2))
    for j in not_in2:
        print("\t%s" %j)
    
    all_diff = set()
    for rxn1 in sbml_1.reactions:
        rxn_id = rxn1.id
        try:
            rxn2 = sbml_2.reactions.get_by_id(rxn_id)
            same_cpd, same_rev =  compare_rxn(rxn1, rxn2)
            if rxn_id not in all_diff:
                if not same_cpd:
                    print("%s use different species" %rxn_id)
                if not same_rev:
                    print("%s use different reversibility" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass
                
    for rxn2 in sbml_2.reactions:
        rxn_id = rxn2.id
        try:
            rxn1 = sbml_1.reactions.get_by_id(rxn_id)
            same_cpd, same_rev =  compare_rxn(rxn1, rxn2)
            if rxn_id not in all_diff:
                if not same_cpd:
                    print("%s use different species" %rxn_id)
                if not same_rev:
                    print("%s use different reversibility" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass

def compare_rxn(rxn1,rxn2):
    """
    compare two cobra reaction object and return (same_cpd, same_rev)
    same_cpd: bool, if true means same compounds consumed and produced
    same_reve: bool, if true means same direction of reaction (reversible or not)

    Parameters
    ----------
    rxn1: cobra.model.reaction
        reaction as cobra object
    rxn2: cobra.model.reaction
        reaction as cobra object

    Returns
    -------
    tuple:
        (same_cpd (bool), same_rev (bool))
    """
    same_cpd, same_rev = False, False
    rxn1_dict,rxn2_dict = {}, {}
    for k,v in list(rxn1.metabolites.items()):
        rxn1_dict[k.id] = v
    for k,v in list(rxn2.metabolites.items()):
        rxn2_dict[k.id] = v
    if rxn1_dict == rxn2_dict: 
        same_cpd = True
    if rxn1.reversibility == rxn2.reversibility:

       same_rev = True

    return (same_cpd, same_rev)

