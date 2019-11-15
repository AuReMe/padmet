# -*- coding: utf-8 -*-
"""
Description:
    For a given padmet file, check and update compartment.
    
    1./ Get all compartment with 1st usage
    
    2./ Remove a compartment with 2nd usage.
    Remove all reactions acting in the given compartment
    
    3./ change compartment id with 3rd usage

"""

def remove_compart(padmet, to_remove, verbose = False):
    """
    Remove all reaction associated to a compound in the compartment to remove.

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    to_remove: str
        compartment id to remove, if many separate compartment id by ','
    verbose: bool
        if True print information

    Returns
    -------
    padmet.classes.PadmetSpec:
        New padmet after removing compartment(s)
    """
    if "," in to_remove:
        to_remove = to_remove.split(",")
    else:
        to_remove = [to_remove]
    for compart in to_remove:
        if verbose:
            print("removing reaction(s) in compartment %s" %compart)
        padmet.delCompart(compart, verbose)
    return padmet

def remplace_compart(padmet, old_compart, new_compart, verbose = False):    
    """
    Remplace compartment 'old_compart' by 'new_compart'.

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    old_comaprt: str
        compartment id to remplace
    new_compart: str
        new compartment id
    verbose: bool
        if True print information

    Returns
    -------
    padmet.classes.PadmetSpec:
        New padmet after remplacing compartment

    """
    if verbose:
        print("Converting compartment %s to %s" %(old_compart, new_compart))
    padmet.change_compart(old_compart, new_compart, verbose)
    return padmet

