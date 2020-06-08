# -*- coding: utf-8 -*-
"""
Description:
    From a given sbml file, create a sbml with only the reactions associated to a gene.

    Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

::

    usage:
        padmet extract_rxn_with_gene_assoc --sbml=FILE --output=FILE [-v]

    options:
        -h --help     Show help.
        --sbml=FILE    path to the sbml file
        --output=FILE    path to the sbml output (with only rxn with genes assoc)
        -v   print info
"""
import docopt
import libsbml
import os

from padmet.utils.sbmlPlugin import parseNotes


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def extract_rxn_with_gene_assoc_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    sbml = args["--sbml"]
    output = args["--output"]
    verbose = args["-v"]

    extract_rxn_with_gene_assoc(sbml, output, verbose)


def extract_rxn_with_gene_assoc(sbml, output, verbose=False):
    """
    From a given sbml document, create a sbml with only the reactions associated to a gene.
    Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

    Parameters
    ----------
    sbml_file: libsbml.document
        sbml document
    output: str
        pathname of the output sbml
    """
    if not os.path.exists(sbml):
        raise FileNotFoundError("No SBML file (--sbml/sbml) accessible at " + sbml)

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML(sbml)
    for i in range(sbml_document.getNumErrors()):
        print(sbml_document.getError(i).getMessage())

    sbml_model = sbml_document.getModel()

    listOfReactions = sbml_model.getListOfReactions()
    
    reactions_to_remove = []
    for reaction in listOfReactions:
        if "GENE_ASSOCIATION" not in list(parseNotes(reaction).keys()):
            reactions_to_remove.append(reaction.getId())
    for rId in reactions_to_remove:
        listOfReactions.remove(rId)
    
    libsbml.writeSBMLToFile(sbml_document, output)

