# -*- coding: utf-8 -*-
"""
Description:
    Allows to visualize a pathway in padmet network.
Color code: 
reactions associated to the pathway, present in the network: lightgreen
reactions associated to the pathway, not present in the network: red
compounds: skyblue

::

    usage:
        padmet visu_path --padmetSpec=FILE/FOLDER --padmetRef=FILE --pathway=ID --output=FILE [--hiding]
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE/FOLDER    pathname to the PADMet file of the network or to a folder containing multiple padmets.
        --padmetRef=FILE    pathname to the PADMet file of the db of reference.
        --pathway=ID    pathway id to visualize, can be multiple pathways separated by a ",".
        --output=FILE    pathname to the output file (extension can be .png or .svg).
        --hiding    hide common compounds.

"""
import docopt
import matplotlib.pylab as plt
import networkx as nx
import os
import seaborn as sns

sns.set_style("white")
sns.set('poster', rc={'figure.figsize':(75,70), 'lines.linewidth': 10})

from collections import OrderedDict
from padmet.classes import PadmetRef, PadmetSpec
from padmet.utils.connection import padmet_to_padmet
from networkx.drawing.nx_agraph import graphviz_layout


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def visu_path_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    padmet_pathname = args["--padmetSpec"]
    padmet_ref_pathname = args["--padmetRef"]
    pathway_ids = args["--pathway"]
    output_file = args["--output"]
    hide_compounds = args["--hiding"]
    visu_path(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file, hide_compounds)


def visu_path(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file, hide_compounds=None):
    """ Extract reactions from pathway and create a comppound/reaction graph.

    Parameters
    ----------
    padmet_pathname: str
        pathname of the padmet file or a folder containing multiple padmet
    padmet_ref_pathname: str
        pathname of the padmetRef file
    pathway_ids: str
        name of the pathway (can be multiple pathways separated by a ',')
    output_file: str
        pathname of the output picture (extension can be .png or .svg)
    hide_compounds: bool
        hide common compounds (like water or proton)
    """
    if os.path.isfile(padmet_pathname):
        padmet = PadmetSpec(padmet_pathname)
    else:
        padmet = padmet_to_padmet.padmet_to_padmet(padmet_pathname)
    padmet_ref = PadmetRef(padmet_ref_pathname)

    pathway_ids = pathway_ids.split(',')
    pwy_all_reactions = []

    if hide_compounds:
        compounds_to_hide = ["WATER", "PROTON", "NAD", "NADH", "ATP", "ADP"]
    else:
        compounds_to_hide = []

    def get_reactions(pathway_id, padmet_ref, pwy_all_reactions):
        all_reactions = [rlt.id_in for rlt in padmet_ref.dicOfRelationOut.get(pathway_id, None)
                        if rlt.type == "is_in_pathway"]
        for reaction_id in all_reactions:
            if reaction_id in padmet_ref.dicOfNode:
                node_reaction = padmet_ref.dicOfNode[reaction_id]
                if node_reaction.type == "pathway":
                    pwy_all_reactions = get_reactions(node_reaction.id, padmet_ref, pwy_all_reactions)
                else:
                    if reaction_id not in pwy_all_reactions:
                        pwy_all_reactions.append(reaction_id)

        return pwy_all_reactions

    for pathway_id in pathway_ids:
        if pathway_id in padmet_ref.dicOfNode:
            tmp_pwy_all_reactions = []
            tmp_pwy_all_reactions = get_reactions(pathway_id, padmet_ref, tmp_pwy_all_reactions)
            pwy_all_reactions.extend(tmp_pwy_all_reactions)
        else:
            print("Pathway " + pathway_id + " not in PadmetRef " + padmet_ref_pathname)

    reactions_in_network = []
    for reaction_id in pwy_all_reactions:
        if reaction_id in padmet.dicOfNode:
            reactions_in_network.append(reaction_id)

    DG=nx.DiGraph()

    custom_node_color = OrderedDict()
    for reaction_id in pwy_all_reactions:
        # Reaction colors
        if reaction_id in reactions_in_network:
            custom_node_color[reaction_id] = "lightgreen"
        else:
            custom_node_color[reaction_id] = "red"

        # Reactants & products for each reaction
        reactants = [rlt.id_out
            for rlt in padmet_ref.dicOfRelationIn.get(reaction_id, None)
                if rlt.type == "consumes"]
        products = [rlt.id_out
            for rlt in padmet_ref.dicOfRelationIn.get(reaction_id, None)
                if rlt.type == "produces"]

        for reac in reactants:
            if reac not in compounds_to_hide:
                if reac not in custom_node_color:
                    custom_node_color[reac] = "skyblue"
                DG.add_edge(reac, reaction_id)
        for prod in products:
            if prod not in compounds_to_hide:
                if prod not in custom_node_color:
                    custom_node_color[prod] = "skyblue"
                DG.add_edge(reaction_id, prod)

    # https://networkx.github.io/documentation/latest/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html
    # apt-get install graphviz graphviz-dev (python-pygraphviz)
    # pip install pygraphviz

    nx.draw_networkx(DG,
                     pos=graphviz_layout(DG, prog='neato'), # Layout from graphviz
                     node_size=1600,
                     arrows=True,
                     font_size=11,      # font-size for labels
                     node_shape='s',    # shape of nodes
                     alpha=0.6,         # node & edge transparency
                     width=1.5,         # line width for edges
                     nodelist=list(custom_node_color.keys()),
                     node_color=[custom_node_color[node] for node in list(custom_node_color.keys())])
    plt.axis('off')

    plt.savefig(output_file, bbox_inches='tight')
    plt.clf()
