#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of padmet.

padmet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
This module contains some functions used for sbml file in addition to libsbml
"""
#pylint: disable=anomalous-backslash-in-string
import re

def parseNotes(element):
    """
    From an SBML element (ex: species or reaction) will return all the section
    note in a dictionary.
    ex:
    <notes>
        <html:body>
            <html:p>BIOCYC: |Alkylphosphonates|</html:p>
            <html:p>CHEBI: 60983</html:p>
        </html:body>
     </notes>
    output: {'BIOCYC': |Alkylphosphonates|,'CHEBI':'60983'}
    value is a list in case diff lines for the same type of info

    @param element: an element from libsbml
    @type element: libsbml.element
    @return: the dictionary of note
    @rtype: dict
    """
    notes = element.getNotesString()
    notesList = notes.splitlines()
    notesDict = {}
    for line in notesList:
        try:
            #line = <html:p>BIOCYC: |Alkylphosphonates|</html:p>
            start = line.index(">")+1
            end = line.index("<", start)
            line = line[start:end]
            #line = BIOCYC: |Alkylphosphonates|
            key, val = line.split(":")
            #line = [BIOCYC,|Alkylphosphonates|]
            key = re.sub(" ", "_", key)
            if len(val) != 0 and val.count(" ") != len(val):
                notesDict[key] = [val]
        except ValueError:
            continue
    return notesDict

def parseGeneAssoc(GeneAssocStr):
    """
    Given a grammar of 'and', 'or' and '(' ')'. Extracts genes ids to a list.
    (geneX and geneY) or geneW' => [geneX,geneY,geneW]
    @param GeneAssocStr: the string containing genes ids
    @type GeneAssocStr: str
    @return: the list of unique ids
    @rtype: list
    """
    #sub '(',')',' ' by ''   sub "and" by "or"
    resultat = re.sub('\(|\)|\s', "", GeneAssocStr).replace("and", "or")
    #create a set by splitting 'or' then convert to list, set for unique genes
    if len(resultat) != 0:
        resultat = list(set(resultat.split("or")))
    else:
        resultat = []
    return resultat

def extractFormula(elementR):
    """
    From an SBML reaction_element will return the formula in a string
    ex: '1.0 FRUCTOSELYSINE_p => 1.0 FRUCTOSELYSINE_c'
    @param elementR: a reaction from libsbml.element
    @type eleemntR: libsbml.element
    @return: the formula
    @rtype: str
    """
    #get direction of reaction
    direction = elementR.getReversible()
    formula = ""
    #generator of reactants
    reactants = [str(reactant.getStoichiometry()) + " " + reactant.getSpecies()
                for reactant in elementR.getListOfReactants()]
    #generator of products
    products = [str(product.getStoichiometry()) + " " + product.getSpecies()
                for product in elementR.getListOfProducts()]

    #formula = ""
    formula = " + ".join(reactants)
    #formula = "1.0 FRUCTOSELYSINE_p + 1.0 Z"
    if direction:
        formula += " <=> "
    else:
        formula += " => "
    formula += " + ".join(products)
    return formula

def convert_to_coded_id(uncoded, _type=None, compart=None):
    """
    convert an id to sbml valid format. First add type of id "R" for reaction
    "M" for compound at the start and the compart at the end.
    _type+"_"+uncoded+"_"+compart
    then replace not allowed char by integer ordinal
    @param uncoded: the original id to code
    @param _type: the type of the id (ex: 'R' or 'M')
    @param _compart: the compartment of the id (ex: 'c' or 'e')
    @type uncoded, _type, _compart: str
    @return: the coded id
    @rtype: str
    """
    #add type and compart
    if _type is not None:
        uncoded = _type+"_"+uncoded
    if compart is not None:
        uncoded += "_"+compart
    #char list that are not allowed in a sbml id
    charlist = ['-', '|', '/', '(', ')', '\'', '=', '#', '*',
                '.', ':', '!', '+', '[', ']', ',', " "]
    for c in charlist:
        #if a banned char in the uncoded id, convert it using the integer ordinal
        uncoded = uncoded.replace(c, "__" + str(ord(c)) + "__")

    return uncoded

def ascii_replace(match):
    """
    recover banned char from the integer ordinal in the reg.match
    """
    return chr(int(match.group(1)))

def convert_from_coded_id(coded):
    """
    convert an id from sbml format to the original id. try to extract the type of
    the id and the compart using strong regular expression
    @param coded: the encoded id
    @type coded: str
    @return: (the uncoded id, type=None, compart=None)
    @rtype: tuple
    """
    #replace DASH from very old sbmls
    coded = coded.replace('_DASH_', '__')
    #an original id starting with int will start with '_' in sbml
    if coded.startswith("_"):
        coded = coded[1:]
    #reg ex to find the ascii used to replace not allowed char
    codepat = re.compile('__(\d+)__')
    #replace ascii by the not allowed char of sbml
    coded = codepat.sub(ascii_replace, coded)

    reg_expr = re.compile('(?P<_type>^[MR]_)(?P<_id>.*)(?P<compart>_.*)')
    search_result = reg_expr.search(coded)
    if search_result is not None:
        compart = search_result.group('compart').replace("_", "")
        _type = search_result.group('_type').replace("_", "")
        uncoded = search_result.group('_id')
    else:
        reg_expr = re.compile('(?P<_type>^[MR]_)(?P<_id>.*)')
        search_result = reg_expr.search(coded)
        if search_result is not None:
            compart = None
            _type = search_result.group('_type').replace("_", "")
            uncoded = search_result.group('_id')
        else:
            reg_expr = re.compile('(?P<_id>.*)(?P<compart>_.*)')
            search_result = reg_expr.search(coded)
            if search_result is not None:
                _type = None
                compart = search_result.group('compart').replace("_", "")
                uncoded = search_result.group('_id')
            else:
                uncoded = coded
                _type = None
                compart = None

    if _type == "R" and compart is not None:
        uncoded += "_" + compart

    return (uncoded, _type, compart)
