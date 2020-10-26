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

    usage:
        padmet padmet_medium --padmetSpec=FILE
        padmet padmet_medium --padmetSpec=FILE -r [--output=FILE] [-v]
        padmet padmet_medium --padmetSpec=FILE --seeds=FILE [--padmetRef=FILE] [--output=FILE] [-v]

    options:
        -h --help     Show help.
        --padmetSpec=FILE    path to the padmet file to update
        --padmetRef=FILE    path to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
        --seeds=FILE    the path to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
        --output=FILE    If not None, pathname to the padmet file updated
        -r    Use to remove all medium from padmet
        -v   print info
"""
import docopt
import os

from padmet.classes import PadmetRef, PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def padmet_medium_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    if args["--seeds"]:
        seeds_file = args["--seeds"]
        if not os.path.exists(seeds_file):
            raise FileNotFoundError("No seeds file (--seeds) accessible at " + seeds_file)

        with open(seeds_file, 'r') as f:
            seeds = [line.split("\t")[0] for line in f.read().splitlines()]
    else:
        seeds = None
    padmet = PadmetSpec(args["--padmetSpec"])
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    output = args["--output"]
    verbose = args["-v"]
    remove = args["-r"]

    if output is None:
        output = args["--padmetSpec"]

    if not remove and not seeds:
        g_m = padmet.get_growth_medium()
        print("List of growth medium:")
        if g_m:
            print(list(g_m))
        else:
            print("[]")
    else:
        manage_medium(padmet, seeds, padmetRef, verbose)
        padmet.generateFile(output)


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

    
    