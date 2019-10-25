# -*- coding: utf-8 -*-
"""
Description:
    From a given sbml file, create a sbml with only the reactions associated to a gene.

    Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

"""
import libsbml
from padmet.utils.sbmlPlugin import parseNotes

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

