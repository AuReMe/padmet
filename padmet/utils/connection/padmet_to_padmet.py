# -*- coding: utf-8 -*-
"""

Description:
    Allows to merge 1-n padmet.
    1./ Update the 'init_padmet' with the 'to_add' padmet(s).
    to_add can be a file or a folder with only padmet files to add.
    
    padmetRef can be use to ensure data uniformization.

"""
from padmet.classes import PadmetSpec
import os

def padmet_to_padmet(to_add, output, padmetRef=None, verbose=False):
    """
    
    """
    if os.path.isdir(to_add):
        path = to_add
        all_files = [i for i in next(os.walk(path))[2] if not i.startswith(".~lock")]
        padmetFiles = [os.path.join(path, i) for i in all_files if i.endswith(".padmet")]
        if len(padmetFiles) == 0:
            print("No padmet found in %s" %path)
            return
    else:
        padmetFiles = to_add.split(";")
    
    padmet_init_file = padmetFiles[0]
    padmet_init = PadmetSpec(padmet_init_file)
    padmetFiles.pop(0)

    for padmet_update_file in padmetFiles:
        if verbose:
            print("Updating %s from %s" %(os.path.basename(padmet_init_file),os.path.basename(padmet_update_file)))
        padmet_update = PadmetSpec(padmet_update_file)
        padmet_init.updateFromPadmet(padmet_update)

    if verbose:
        print("Generated file: %s" %output)
    padmet_init.generateFile(output)
        
