# -*- coding: utf-8 -*-
"""
Description:
    For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
    for each compounds to maintain consistency of the network for flux analysis.
    For each compounds create:

        An exchange reaction: this reaction consumes the compound in the
        compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular

        A transport reaction: this reaction consumes the compound in the compartment
        'e' extracellular' and produces the compound in the compartment 'c' cytosol
        ex: for seed 'cpd-a'

    1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.

    2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) <=> 1 cpd-a (e)

    3/ create transport reaction: TransportSeed_cpd-a_e: 1 cpd-a (e) => 1 cpd-a (c)

    4/ create a new file if output not None, or overwrite padmetSpec

::
"""
def manage_medium(padmet, new_growth_medium=None, padmetRef=None, verbose=False):
    """
    Manage medium of a padmet. If new_growth_medium give, use this list of compound
    to define the new medium and create transport and exchange reactions.
    if padmetRef given, use the information from padmetRef to create the missing compound.
    If no new_growth_medium given: remove the current medium in the padmet.

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to update
    new_growth_medium: list
        list of compound id representing the medium
    padmetRef: padmet.classes.PadmetRef
        padmet containing the database of reference
    verbose: bool
        if True print information

    Returns
    -------
    padmet.classes.PadmetSpec:
        New padmet after updating medium
    """
    if new_growth_medium:
        padmet.set_growth_medium(new_growth_medium=new_growth_medium, padmetRef=padmetRef, verbose=verbose)
    else:
        padmet.remove_growth_medium(verbose=verbose)
    return padmet

    
    