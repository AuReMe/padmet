# -*- coding: utf-8 -*-
"""
Description:
    From the padmetRef of MetaCyc creates the MetaCyc ontology.
    At this moment, all the element of the tree begins with a ont_ and all the '+' or '-' are removed.
    This is a limitation from lxml tag.

::

    usage:
        padmet get_metacyc_to_ontology -p=FILE -o=FILE

    options:
        -h --help     Show help.
        -p=FILE    path of the padmet file of MetaCyc
        -o=FILE   pathname of the XML output file
"""
import docopt
from lxml import etree
from padmet.classes.padmetRef import PadmetRef


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def get_metacyc_ontology_cli(command_args):
    #parsing args
    args = docopt.docopt(__doc__, argv=command_args)
    padmetRef_file = args["-p"]
    output_file = args["-o"]
    metacyc_to_ontology(padmetRef_file, output_file, ontology_root='FRAMES')


def metacyc_to_ontology(padmetRef_file, output_file, ontology_root='FRAMES'):
    """
    Extract teh ontology of MetaCyc from the padmetRef.

    Parameters
    ----------
    padmetRef_file: str
        path to padmetRef file
    output_file: str
        pathname of the output sbml
    ontology_root: str
        name of the roots to use to create the tree (FRAMES, Generalized-Reactions, Compounds, ...)
    """
    padmetref = PadmetRef(padmetRef_file)

    class_nodes = [node for node in padmetref.dicOfNode.values() if node.type == "class"]

    known_parents = {}

    # Create the root of the xml tree.
    frames = etree.Element("element_" + str(len(known_parents)), name=ontology_root)

    known_parents[ontology_root] = frames

    # Extract the parent of each objects from classes.
    child_parents = {}
    for class_node in class_nodes:
        if class_node.id != 'FRAMES':
            parent_classes = [rlt.id_out for rlt in padmetref.dicOfRelationIn[class_node.id] if rlt.type == 'is_a_class']
            child_parents[class_node.id] = parent_classes

    def get_child(parent_id, known_parents, child_parents):
        # Search all the child of the parent_id
        for child_parent in child_parents:
            for parent_class in child_parents[child_parent]:
                if parent_class == parent_id:
                    if parent_class in known_parents:
                        et_subelement = etree.SubElement(known_parents[parent_class], "element_" + str(len(known_parents)), name=child_parent)
                        known_parents[child_parent] = et_subelement
                        get_child(child_parent, known_parents, child_parents)


    get_child(ontology_root, known_parents, child_parents)

    tree = etree.ElementTree(frames)
    tree.write(output_file, pretty_print=True)


def add_element_to_tree(element_ids, padmet_instance, ontology_elements, element_type):
    count = 0
    for element_id in element_ids:
        if element_id in padmet_instance.dicOfRelationIn:
            if element_type == 'reaction':
                element_classes = [rlt.id_out for rlt in padmet_instance.dicOfRelationIn[element_id] if rlt.type == "is_in_pathway"]
            else:
                element_classes = [rlt.id_out for rlt in padmet_instance.dicOfRelationIn[element_id] if rlt.type == "is_a_class"]
            for element_class in element_classes:
                if element_class in ontology_elements:
                    for subclass in ontology_elements[element_class]:
                        etree_sublement = etree.SubElement(subclass, element_type + '_element_' + str(count), name=element_id)
                        if element_id not in ontology_elements:
                            ontology_elements[element_id] = [etree_sublement]
                        else:
                            ontology_elements[element_id].append(etree_sublement)
                        count += 1

    return ontology_elements


def extract_element_ontology(metacyc_ontology_file, padmetRef_file, output_file):
    onttree = etree.parse(metacyc_ontology_file)
    ontology_elements = {}
    for element in onttree.iter():
        if element.attrib['name'] not in ontology_elements:
            ontology_elements[element.attrib['name']] = [element]
        else:
            ontology_elements[element.attrib['name']].append(element)

    padmetref = PadmetRef(padmetRef_file)

    compound_ids = [node.id for node in padmetref.dicOfNode.values() if node.type == "compound"]
    ontology_elements = add_element_to_tree(compound_ids, padmetref, ontology_elements, 'compound')

    pathway_ids = [node.id for node in padmetref.dicOfNode.values() if node.type == "pathway"]
    ontology_elements = add_element_to_tree(pathway_ids, padmetref, ontology_elements, 'pathway')

    reaction_ids = [node.id for node in padmetref.dicOfNode.values() if node.type == "reaction"]
    ontology_elements = add_element_to_tree(reaction_ids, padmetref, ontology_elements, 'reaction')

    tree = etree.ElementTree(onttree.getroot())
    tree.write(output_file, pretty_print=True)


def ontology_to_newick(metacyc_ontology_file, newick_output_file):
    onttree = etree.parse(metacyc_ontology_file)

    def child_to_newick(node):
        childrens = []
        if node.getparent() is None:
            childrens.append(node.attrib['name'])
        for children in node.getchildren():
            if len(children.getchildren()) > 0:
                subchilds = []
                childs = child_to_newick(children)
                subchilds.append("'"+children.attrib['name']+"'")
                subchilds.append(childs)
                childrens.append('(' + ','.join(subchilds) + ')')
            else:
                childrens.append("'"+children.attrib['name']+"'")
        return '(' + ','.join(childrens) + ')'

    with open(newick_output_file, 'w') as newick:
        newick.write(child_to_newick(onttree.getroot())+';')
