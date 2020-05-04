# -*- coding: utf-8 -*-
"""
Description:
    Create reports of a padmet file. 

    all_pathways.tsv: header = ["dbRef_id", "Common name", "Number of reaction found",
    "Total number of reaction", "Ratio (Reaction found / Total)"]

    all_reactions.tsv: header = ["dbRef_id", "Common name", "formula (with id)",
    "formula (with common name)", "in pathways", "associated genes"]
    
    all_metabolites.tsv: header = ["dbRef_id", "Common name", "Produced (p), Consumed (c), Both (cp)"]

::
    
    usage:
        padmet report_network --padmetSpec=FILE --output_dir=dir [--padmetRef=FILE] [-v]
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE    pathname of the padmet file.
        --padmetRef=FILE    pathname of the padmet file used as database
        --output_dir=dir    directory for the results.
        -v   print info.
"""
import docopt

from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def report_network_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    padmetRef_file = args["--padmetRef"]
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    output_dir = args["--output_dir"]
    verbose = args["-v"]

    padmetSpec.network_report(output_dir, padmetRef_file, verbose)  
