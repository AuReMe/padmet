# -*- coding: utf-8 -*-
"""

Description:
    Allows to merge 1-n padmet.
    1./ Update the 'init_padmet' with the 'to_add' padmet(s).
    to_add can be a file or a folder with only padmet files to add.

::

    usage:
        padmet padmet_to_padmet --to_add=FILE/DIR --output=FILE  [-v]

    options:
        -h --help     Show help.
        --to_add=FILE/DIR    path to the padmet file to add (sep: ,) or path to folder of padmet files.
        --output=FILE   path to the new padmet file
        -v   print info
"""
import docopt
import os

from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def padmet_to_padmet_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    output = args["--output"]
    verbose = args["-v"]
    to_add = args["--to_add"]
    padmet_to_padmet(to_add, output, verbose)


def padmet_to_padmet(to_add, output=None, verbose=False):
    """
    Create a padmet by merging multiple other padmet files.

    Parameters
    ----------
    to_add: dir or str
        padmet directory or string with multiple padmet paths separated by ','
    output:
        path to the output file
    verbose: bool
        verbose level of script
    Returns
    -------
    padmet_init: PadmetSpec
        padmet created from emrging of the other padmet
    """
    if os.path.isdir(to_add):
        path = to_add
        all_files = [i for i in next(os.walk(path))[2] if not i.startswith(".~lock")]
        padmetFiles = [os.path.join(path, i) for i in all_files if i.endswith(".padmet")]
        if len(padmetFiles) == 0:
            print("No padmet found in %s" %path)
            return
    else:
        padmetFiles = to_add.split(",")
    
    padmet_init_file = padmetFiles[0]
    padmet_init = PadmetSpec(padmet_init_file)
    padmetFiles.pop(0)

    for padmet_update_file in padmetFiles:
        if verbose:
            print("Updating %s from %s" %(os.path.basename(padmet_init_file),os.path.basename(padmet_update_file)))
        padmet_update = PadmetSpec(padmet_update_file)
        padmet_init.updateFromPadmet(padmet_update)

    if output:
        if verbose:
            print("Generated file: %s" %output)
        padmet_init.generateFile(output)

    return padmet_init