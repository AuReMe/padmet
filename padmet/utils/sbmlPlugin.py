#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=anomalous-backslash-in-string
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

    Parameters
    ----------
    element: libsbml.element
        an element from libsbml

    Returns
    -------
    dict:
        the dictionary of note
    """
    notes = element.getNotesString()
    notes_list = notes.splitlines()
    notes_dict = {}
    for line in notes_list:
        try:
            # line = <html:p>BIOCYC: |Alkylphosphonates|</html:p>
            start = line.index(">") + 1
            end = line.index("<", start)
            line = line[start:end]
            # line = BIOCYC: |Alkylphosphonates|
            key, val = line.split(":")
            # line = [BIOCYC,|Alkylphosphonates|]
            key = re.sub(" ", "_", key)
            if len(val) != 0 and val.count(" ") != len(val):
                notes_dict[key] = [val]
        except ValueError:
            continue
    return notes_dict


def parseGeneAssoc(GeneAssocStr):
    """
    Given a grammar of 'and', 'or' and '(' ')'. Extracts genes ids to a list.
    (geneX and geneY) or geneW' => [geneX,geneY,geneW]
    Parameters
    ----------
    GeneAssocStr: str
        the string containing genes ids

    Returns
    -------
    list:
        the list of unique ids
    """
    # remplace ' and ' or ' or ' by a tag '_FORSPLIT_'
    GeneAssocStr_tmp = re.sub(r" and | or ", "_FORSPLIT_", GeneAssocStr)
    # remove '(' or ')' or ' '
    resultat = re.sub(r"\(|\)|\s", "", GeneAssocStr_tmp)
    # create a set by splitting '_FORSPLIT_' then convert to list, set for unique genes
    if len(resultat) != 0:
        resultat = list(set(resultat.split("_FORSPLIT_")))
    else:
        resultat = []
    return resultat


def extractFormula(elementR):
    """
    From an SBML reaction_element will return the formula in a string
    ex: '1.0 FRUCTOSELYSINE_p => 1.0 FRUCTOSELYSINE_c'
    Parameters
    ----------
    elementR: libsbml.element
        a reaction from libsbml.element

    Returns
    -------
    str:
        the formula
    """
    # get direction of reaction
    direction = elementR.getReversible()
    formula = ""
    # generator of reactants
    reactants = [
        str(reactant.getStoichiometry()) + " " + reactant.getSpecies()
        for reactant in elementR.getListOfReactants()
    ]
    # generator of products
    products = [
        str(product.getStoichiometry()) + " " + product.getSpecies()
        for product in elementR.getListOfProducts()
    ]

    # formula = ""
    formula = " + ".join(reactants)
    # formula = "1.0 FRUCTOSELYSINE_p + 1.0 Z"
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
    Parameters
    ----------
    uncoded: str
        the original id to code
    _type: str
        the type of the id (ex: 'R' or 'M')
    compart: str
        the compartment of the id (ex: 'c' or 'e')

    Returns
    -------
    str:
        the coded id
    """
    # add type and compart
    if _type is not None:
        uncoded = _type + "_" + uncoded
    if compart is not None:
        uncoded += "_" + compart
    # char list that are not allowed in a sbml id
    charlist = [
        "-",
        "|",
        "/",
        "(",
        ")",
        "'",
        "=",
        "#",
        "*",
        ".",
        ":",
        "!",
        "+",
        "[",
        "]",
        ",",
        " ",
    ]
    for char in charlist:
        # if a banned char in the uncoded id, convert it using the integer ordinal
        uncoded = uncoded.replace(char, "__" + str(ord(char)) + "__")

    return uncoded


def ascii_replace(match):
    """
    recover banned char from the integer ordinal in the reg.match
    """
    return chr(int(match.group(1)))


def convert_from_coded_id(
    coded, pattern="__", compart_in_id=False, reaction_tag="R", species_tag="M"
):
    """
    convert an id from sbml format to the original id. try to extract the type of
    the id and the compart using strong regular expression
    ex: M_METABOLITE__45__12_c => ('METABOLITE-12', 'M', 'c')

    Parameters
    ----------
    coded: str
        the encoded id
    pattern: str
        pattern used to delimit interger ordinal
    compart_in_id: bool
        if true: the last _* is not mean to be the compart is part of the id
    reaciton_tag: str
        First letter used to tag a reaction
    species_tag: str
        First letter used to tag a species

    Returns
    -------
    str:
        the uncoded id
    str, None:
        type of ID (ex: 'M' or 'R')
    str, None:
        compart of the id
    """
    # replace DASH from very old sbmls
    coded = coded.replace("_DASH_", "__")
    # an original id starting with int will start with '_' in sbml
    if coded.startswith("_"):
        coded = coded[1:]
    # reg ex to find the ascii used to replace not allowed char
    ascii_pattern = r"{0}(\d+){0}".format(pattern)
    codepat = re.compile(ascii_pattern)
    # replace ascii by the not allowed char of sbml
    coded = codepat.sub(ascii_replace, coded)

    # Regular expression for coded ID from:
    # https://github.com/SBRG/bigg_models/wiki/BiGG-Models-ID-Specification-and-Guidelines
    str_reg = r"(?P<_type>^[{0}{1}]_)(?P<_id>.*)(?P<compart>_[a-z][a-z0-9]?$)".format(
        species_tag, reaction_tag
    )
    reg_expr = re.compile(str_reg)
    search_result = reg_expr.search(coded)
    if search_result is not None:
        compart = search_result.group("compart").replace("_", "")
        _type = search_result.group("_type").replace("_", "")
        uncoded = search_result.group("_id")
    else:
        str_reg = r"(?P<_type>^[{0}{1}]_)(?P<_id>.*)".format(species_tag, reaction_tag)
        reg_expr = re.compile(str_reg)
        search_result = reg_expr.search(coded)
        if search_result is not None:
            compart = None
            _type = search_result.group("_type").replace("_", "")
            uncoded = search_result.group("_id")
        else:
            reg_expr = re.compile(r"(?P<_id>.*)(?P<compart>_.*)")
            search_result = reg_expr.search(coded)
            if search_result is not None:
                _type = None
                compart = search_result.group("compart").replace("_", "")
                uncoded = search_result.group("_id")
            else:
                uncoded = coded
                _type = None
                compart = None

    if compart and compart_in_id:
        uncoded += "_" + compart

    return (uncoded, _type, compart)


def get_all_decoded_version(element_id, _type):
    """
    Use convert_from_coded function to convert a element_id (reaction or species)
    _type use define if element is a 'reaction' or un 'species'.
    Try different decoding combination based on old and new sbml id encoding.

    Parameters
    ----------
    element_id: str
        the encoded id
    _type: str
        _type is 'reaction' or 'species'

    Returns
    -------
    list:
        list of encoded id
    """
    all_element_id_decoded = list()
    # 1st attemp: decoded id with classic encoding convention. cf sbmlplugin.convert_from_coded_id
    all_element_id_decoded.append(convert_from_coded_id(element_id)[0])
    if _type == "species":
        # 2st attemp: decoded id with classic non-conventionnal encoding. cf sbmlplugin.convert_from_coded_id
        if (
            convert_from_coded_id(element_id, compart_in_id=True)[0]
            not in all_element_id_decoded
        ):
            all_element_id_decoded.append(
                convert_from_coded_id(element_id, compart_in_id=True)[0]
            )
        if (
            convert_from_coded_id(element_id, pattern="_", species_tag="S")[0]
            not in all_element_id_decoded
        ):
            all_element_id_decoded.append(
                convert_from_coded_id(element_id, pattern="_", species_tag="S")[0]
            )
    elif _type == "reaction":
        # 1st attemp: decoded id with classic encoding convention. cf sbmlplugin.convert_from_coded_id
        all_element_id_decoded.append(convert_from_coded_id(element_id)[0])
        # 2st attemp: decoded id with classic non-conventionnal encoding. cf sbmlplugin.convert_from_coded_id
        if (
            convert_from_coded_id(element_id, compart_in_id=True)[0]
            not in all_element_id_decoded
        ):
            all_element_id_decoded.append(
                convert_from_coded_id(element_id, compart_in_id=True)[0]
            )
    return all_element_id_decoded

