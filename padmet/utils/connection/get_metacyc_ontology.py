# -*- coding: utf-8 -*-
"""
Description:
    From the padmetRef of MetaCyc creates the MetaCyc ontology.
    At this moment, all the element of the tree begins with a ont_ and all the '+' or '-' are removed.
    This is a limitation from lxml tag.
"""
from lxml import etree
from padmet.classes.padmetRef import PadmetRef

def metacyc_to_ontology(padmetRef_file, output_file):
    """
    Extract teh ontology of MetaCyc from the padmetRef.

    Parameters
    ----------
    padmetRef_file: str
        path to padmetRef file
    output_file: str
        pathname of the output sbml
    """
    padmetref = PadmetRef(padmetRef_file)

    class_nodes = [node for node in padmetref.dicOfNode.values() if node.type == "class"]

    known_parents = {}

    # FRAMES can be replaced wiht Generalized-Reactions (for reactions/pathways) or Compounds
    frames = etree.Element("FRAMES")

    known_parents['FRAMES'] = frames

    child_parents = {}

    for class_node in class_nodes:
        if class_node.id != 'FRAMES':
            parent_classes = [rlt.id_out for rlt in padmetref.dicOfRelationIn[class_node.id] if rlt.type == 'is_a_class']
            child_parents[class_node.id] = parent_classes

    def get_child(parent_id, known_parents, child_parents):
        for child_parent in child_parents:
            for parent_class in child_parents[child_parent]:
                if parent_class == parent_id:
                    if parent_class in known_parents:
                        et_subelement = etree.SubElement(known_parents[parent_class], 'ont_'+child_parent.replace('+','').replace(':',''))
                        known_parents[child_parent] = et_subelement
                        get_child(child_parent, known_parents, child_parents)


    get_child('FRAMES', known_parents, child_parents)

    tree = etree.ElementTree(frames)
    tree.write(output_file, pretty_print=True)