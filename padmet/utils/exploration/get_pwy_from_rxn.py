# -*- coding: utf-8 -*-
"""
Description:
    From a file containing a list of reaction, return the pathways where these reactions 
    are involved.
    ex: if rxn-a in pwy-x => return, pwy-x; all rxn ids in pwy-x; all rxn ids in pwy-x FROM the list; ratio

"""

def extract_pwys(padmet, reactions):
    """
    #extract from padmet pathways containing 1-n reactions from a set of reactions 'reactions'
    Return a dict of data.
    dict, k=pathway_id, v=dict: k in [total_rxn, rxn_from_list, ratio
    ex: {pwy-x:{'total_rxn':[a,b,c], rxn_from_list:[a], ratio:1/3}}


    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    reactions: set
        set of reactions to match with pathways
        
    Returns
    -------
    dict:
        dict, k=pathway_id, v=dict: k in [total_rxn, rxn_from_list, ratio
        ex: {pwy-x:{'total_rxn':[a,b,c], rxn_from_list:[a], ratio:1/3}}
    """
    dict_pwy = {}
    pwys_id = [node.id for node in padmet.dicOfNode.values() if node.type == "pathway"]
    for pwy_id in pwys_id:
        rxn_in_pwy = set([rlt.id_in for rlt in padmet.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway" and padmet.dicOfNode[rlt.id_in].type == "reaction"])
        rxn_from_list = reactions.intersection(rxn_in_pwy)
        if rxn_from_list:
            ratio = round(float(len(rxn_from_list))/float(len(rxn_in_pwy)),2)
            dict_pwy[pwy_id] = {'total_rxn': rxn_in_pwy, 'rxn_from_list': rxn_from_list, 'ratio': ratio}
    return dict_pwy
    
def dict_pwys_to_file(dict_pwy, output):
    """
    Create csv file from dict_pwy. dict_pwy is obtained with extract_pwys()

    Parameters
    ----------
    dict_pwy: dict
        dict, k=pathway_id, v=dict: k in [total_rxn, rxn_from_list, ratio
        ex: {pwy-x:{'total_rxn':[a,b,c], rxn_from_list:[a], ratio:1/3}}
    output: str
        path to output file

    """
    with open(output, 'w') as f:
        header = ["pathway_id","total_rxn","rxn_from_list","ratio"]
        f.write("\t".join(header)+"\n")
        for pwy_id, pwy_data in dict_pwy.items():
            total_rxn = pwy_data['total_rxn']
            rxn_from_list = pwy_data['rxn_from_list']
            ratio = pwy_data['ratio']

            line = [pwy_id, ";".join(total_rxn), ";".join(rxn_from_list),str(ratio)]
            f.write("\t".join(line)+"\n")
