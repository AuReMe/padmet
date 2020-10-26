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
        padmet visu_path --padmetSpec=FILE/FOLDER --padmetRef=FILE --pathway=ID --output=FILE [--hide-currency] [--level=STR]
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE/FOLDER    pathname to the PADMet file of the network or to a folder containing multiple padmets.
        --padmetRef=FILE    pathname to the PADMet file of the db of reference.
        --pathway=ID    pathway id to visualize, can be multiple pathways separated by a ",".
        --output=FILE    pathname to the output file (extension can be .png or .svg).
        --hide-currency    hide currency metabolites.
        --level=STR    level of precision for the visualization (compound or pathway). By default visualization uses "compound".
"""
import docopt
import matplotlib.pylab as plt
import networkx as nx
import os
import seaborn as sns
import sys

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
    hide_currency_metabolites = args["--hide-currency"]
    visualization_level = args["--level"]
    if visualization_level is None:
        visualization_level = "compound"
    if visualization_level == "compound":
        visu_path_compounds(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file, hide_currency_metabolites)
    if visualization_level == "pathway":
        visu_path_pathways(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file)


def visu_path_compounds(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file, hide_currency_metabolites=None):
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
    hide_currency_metabolites: bool
        hide currency metabolites
    """
    if os.path.isfile(padmet_pathname):
        padmet = PadmetSpec(padmet_pathname)
    else:
        padmet = padmet_to_padmet.padmet_to_padmet(padmet_pathname)
    padmet_ref = PadmetRef(padmet_ref_pathname)

    pathway_ids = pathway_ids.split(',')
    pwy_all_reactions = []

    if hide_currency_metabolites:
        compounds_to_hide = ["PROTON", "WATER", "OXYGEN-MOLECULE", "NADP", "NADPH", "ATP",
                            "PPI", "CARBON-DIOXIDE", "Pi", "ADP", "CO-A", "UDP", "NAD",
                            "NADH", "AMP", "AMMONIA", "HYDROGEN-PEROXIDE", "Acceptor",
                            "Donor-H2", "3-5-ADP", "GDP", "CARBON-MONOXIDE", "GTP", "FAD"]
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
                if 'REVERSIBLE' in padmet_ref.dicOfNode[reaction_id].misc['DIRECTION']:
                    DG.add_edge(reaction_id, reac)
        for prod in products:
            if prod not in compounds_to_hide:
                if prod not in custom_node_color:
                    custom_node_color[prod] = "skyblue"
                DG.add_edge(reaction_id, prod)
                if 'REVERSIBLE' in padmet_ref.dicOfNode[reaction_id].misc['DIRECTION']:
                    DG.add_edge(prod, reaction_id)

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


def visu_path_pathways(padmet_pathname, padmet_ref_pathname, pathway_ids, output_file):
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

    # Check if the padmets and padmetref contain the INPUT-COMPOUNDS and OUTPUT-COMPOUNDS in pathway node.misc needed for this analysis.
    padmetref_input_compounds_in_pwys = [1 for node_pathway in padmet_ref.dicOfNode
                                            if padmet_ref.dicOfNode[node_pathway].type == 'pathway' and 'INPUT-COMPOUNDS' in padmet_ref.dicOfNode[node_pathway].misc]
    padmetref_output_compounds_in_pwys = [1 for node_pathway in padmet_ref.dicOfNode
                                            if padmet_ref.dicOfNode[node_pathway].type == 'pathway' and 'OUTPUT-COMPOUNDS' in padmet_ref.dicOfNode[node_pathway].misc]
    if sum(padmetref_input_compounds_in_pwys) == 0 or sum(padmetref_output_compounds_in_pwys) == 0:
        sys.exit("The padmetref " + padmet_ref_pathname + " does not contain INPUT-COMPOUNDS and OUTPUT-COMPOUNDS in the pathway node, can't produce the pathway visualization.")

    padmet_input_compounds_in_pwys = [1 for node_pathway in padmet.dicOfNode
                                        if padmet.dicOfNode[node_pathway].type == 'pathway' and 'INPUT-COMPOUNDS' in padmet.dicOfNode[node_pathway].misc]
    padmet_output_compounds_in_pwys = [1 for node_pathway in padmet.dicOfNode
                                        if padmet.dicOfNode[node_pathway].type == 'pathway' and 'OUTPUT-COMPOUNDS' in padmet.dicOfNode[node_pathway].misc]
    if sum(padmet_input_compounds_in_pwys) == 0 or sum(padmet_output_compounds_in_pwys) == 0:
        sys.exit("The padmet " + padmet_pathname + " does not contain INPUT-COMPOUNDS and OUTPUT-COMPOUNDS in the pathway node, can't produce the pathway visualization.")

    # Extract pathway from superpathways.
    pathway_ids = pathway_ids.split(',')
    all_pathways = []

    def get_pathways(pathway_id, padmet_ref, pwy_all_reactions):
        all_reactions_pathways = [rlt.id_in for rlt in padmet_ref.dicOfRelationOut.get(pathway_id, None)
                        if rlt.type == "is_in_pathway"]
        for reaction_pathway_id in all_reactions_pathways:
            if reaction_pathway_id in padmet_ref.dicOfNode:
                node_reaction = padmet_ref.dicOfNode[reaction_pathway_id]
                if node_reaction.type == "pathway":
                    pwy_all_reactions.append(reaction_pathway_id)
                    pwy_all_reactions = get_pathways(node_reaction.id, padmet_ref, pwy_all_reactions)

        return pwy_all_reactions

    for pathway_id in pathway_ids:
        if pathway_id in padmet_ref.dicOfNode:
            tmp_pwy_all_pathways = []
            tmp_pwy_all_pathways = get_pathways(pathway_id, padmet_ref, tmp_pwy_all_pathways)
            all_pathways.extend(tmp_pwy_all_pathways)
        else:
            print("Pathway " + pathway_id + " not in PadmetRef " + padmet_ref_pathname)

    # Find pathway in the padmet file.
    pathways_in_network = []
    for pathway_id in all_pathways:
        if pathway_id in padmet.dicOfNode:
            pathways_in_network.append(pathway_id)

    # Create the graph.
    DG=nx.DiGraph()
    custom_node_color = OrderedDict()

    for pwy in all_pathways:
        node_pathway = padmet_ref.dicOfNode[pwy]
        if pwy in pathways_in_network:
            custom_node_color[pwy] = "lightgreen"
        else:
            custom_node_color[pwy] = "red"
        if 'INPUT-COMPOUNDS' in node_pathway.misc and 'OUTPUT-COMPOUNDS' in node_pathway.misc:
            for reactant in node_pathway.misc['INPUT-COMPOUNDS'][0].split(','):
                if reactant not in custom_node_color:
                    custom_node_color[reactant] = "skyblue"
                DG.add_edge(reactant, pwy)
            for product in node_pathway.misc['OUTPUT-COMPOUNDS'][0].split(','):
                if product not in custom_node_color:
                    custom_node_color[product] = "skyblue"
                DG.add_edge(pwy, product)

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
