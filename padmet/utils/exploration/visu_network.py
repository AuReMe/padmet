# -*- coding: utf-8 -*-
"""
Description:
    Allows to visualize a metabolic network on a compounds perspectives

::

    usage:
        padmet visu_network -i=FILE -o=FILE [--html=FILE] [--level=STR] [--hide-currency]
    
    options:
        -h --help     Show help.
        -i=FILE    pathname to the input file (either PADMet or SBML).
        -o=FILE    pathname to the output file (picture of metabolic network).
        --html=FILE    pathname to the output file (interactive hmtl of metabolic network).
        --level=STR    level of precision for the visualization (compound, reaction or pathway). By default visualization uses "compound".
        --hide-currency    hide currency metabolites.
"""
import docopt
import logging
import sys

try:
    import igraph
except ImportError:
    raise ImportError('Requires igraph with cairocffi, try:\npip install python-igraph cairocffi')

from cobra.io.sbml import read_sbml_model
from padmet.classes import PadmetRef, PadmetSpec
from padmet.utils.sbmlPlugin import convert_from_coded_id

logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def visu_network_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    
    metabolic_network_file = args["-i"]
    output_file = args["-o"]
    html_output_file = args['--html']
    visualization_level = args['--level']
    hide_currency_metabolites = args["--hide-currency"]
    if visualization_level is None:
        visualization_level = "compound"
    create_graph(metabolic_network_file, output_file, visualization_level, hide_currency_metabolites)
    if html_output_file:
        create_html_graph(metabolic_network_file, html_output_file, visualization_level)


def parse_compounds_padmet(padmet_file, hide_metabolites):
    """ Parse padmets files to extract compounds to create edges and nodes for igraph.

    Parameters
    ----------
    padmet_file: str
        pathname of the padmet file
    hide_metabolites: list
        list of metabolites to hide
    Returns
    -------
    edges: list
        edges between two compounds (symbolizing the reaction)
    edges_label: list
        for each edge the name of the reaction
    weights: list
        the weight associated to each edge
    nodes: list
        a compound
    nodes_label: list
        for each node the name of the compound
    """
    padmetSpec = PadmetSpec(padmet_file)

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]

    edges = []
    edges_label = []
    weights = []

    nodes = {}
    nodes_label = []
    for rxn in all_rxns:
        ins = []
        outs = []
        for rlt in padmetSpec.dicOfRelationIn[rxn.id]:
            if rlt.type == "consumes":
                ins.append(rlt.id_out)
            if rlt.type == "produces":
                outs.append(rlt.id_out)
        for compound_in in ins:
            if compound_in not in hide_metabolites:
                for compound_out in outs:
                    if compound_out not in hide_metabolites:
                        if compound_in not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(compound_in)
                            nodes[compound_in] = new_cpd_id
                        if compound_out not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(compound_out)
                            nodes[compound_out] = new_cpd_id
                        edges.append((nodes[compound_in], nodes[compound_out]))
                        weights.append(1)
                        edges_label.append(rxn.id)
                        if 'REVERSIBLE' in rxn.misc['DIRECTION']:
                            edges.append((nodes[compound_out], nodes[compound_in]))
                            weights.append(1)
                            edges_label.append(rxn.id)

    return edges, edges_label, weights, nodes, nodes_label


def parse_compounds_sbml(sbml_file, hide_metabolites):
    """ Parse sbml files to extract compounds to create edges and nodes for igraph.

    Parameters
    ----------
    sbml_file: str
        pathname of the sbml file
    hide_metabolites: list
        list of metabolites to hide
    Returns
    -------
    edges: list
        edges between two compounds (symbolizing the reaction)
    edges_label: list
        for each edge the name of the reaction
    weights: list
        the weight associated to each edge
    nodes: list
        a compound
    nodes_label: list
        for each node the name of the compound
    """
    sbml_model = read_sbml_model(sbml_file)

    edges = []
    edges_label = []
    weights = []

    nodes = {}
    nodes_label = []

    for reaction in sbml_model.reactions:
        for reactant in reaction.reactants:
            reactant = convert_from_coded_id(reactant.id)[0]
            if reactant not in hide_metabolites:
                for product in reaction.products:
                    product = convert_from_coded_id(product.id)[0]
                    if product not in hide_metabolites:
                        if reactant not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(reactant)
                            nodes[reactant] = new_cpd_id
                        if product not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(product)
                            nodes[product] = new_cpd_id
                        edges.append((nodes[reactant], nodes[product]))
                        weights.append(1)
                        edges_label.append(reaction.id)
                        if reaction.reversibility == True:
                            edges.append((nodes[product], nodes[reactant]))
                            weights.append(1)
                            edges_label.append(reaction.id)

    return edges, edges_label, weights, nodes, nodes_label


def parse_reactions_padmet(padmet_file):
    """ Parse padmets files to extract reactions to create edges and nodes for igraph.

    Parameters
    ----------
    padmet_file: str
        pathname of the padmet file
    Returns
    -------
    edges: list
        edges between two reactions
    edges_label: list
        for each edge the name of the reaction
    weights: list
        the weight associated to each edge
    nodes: list
        a compound
    nodes_label: list
        for each node the name of the compound
    """
    padmetSpec = PadmetSpec(padmet_file)

    edges = []
    edges_label = []
    weights = []

    nodes = {}
    nodes_label = []

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]

    for rxn in all_rxns:
        for rlt in padmetSpec.dicOfRelationIn[rxn.id]:
            if rlt.type == "produces":
                for sec_rlt in padmetSpec.dicOfRelationOut[rlt.id_out]:
                    if sec_rlt.type == "consumes":
                        sec_rxn_id = sec_rlt.id_in
                        if sec_rxn_id not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(sec_rxn_id)
                            nodes[sec_rxn_id] = new_cpd_id
                        if rxn.id not in nodes:
                            new_cpd_id = len(nodes_label)
                            nodes_label.append(rxn.id)
                            nodes[rxn.id] = new_cpd_id
                        if (nodes[rxn.id], nodes[sec_rxn_id]) not in edges:
                            edges.append((nodes[rxn.id], nodes[sec_rxn_id]))
                            weights.append(1)
                            edges_label.append(rxn.id)

    return edges, edges_label, weights, nodes, nodes_label


def parse_pathways_padmet(padmet_file):
    """ Parse padmets files to extract pathway inputs and ouputs to create edges and nodes for igraph.

    Parameters
    ----------
    padmet_file: str
        pathname of the padmet file
    Returns
    -------
    edges: list
        edges between two compounds (symbolizing the pathway)
    edges_label: list
        for each edge the name of the pathway
    weights: list
        the weight associated to each edge
    nodes: list
        a compound
    nodes_label: list
        for each node the name of the compound
    """
    padmetSpec = PadmetSpec(padmet_file)

    # Check if the padmets and padmetref contain the INPUT-COMPOUNDS and OUTPUT-COMPOUNDS in pathway node.misc needed for this analysis.
    padmet_input_compounds_in_pwys = [1 for node_pathway in padmetSpec.dicOfNode
                                            if padmetSpec.dicOfNode[node_pathway].type == 'pathway' and 'INPUT-COMPOUNDS' in padmetSpec.dicOfNode[node_pathway].misc]
    padmet_output_compounds_in_pwys = [1 for node_pathway in padmetSpec.dicOfNode
                                            if padmetSpec.dicOfNode[node_pathway].type == 'pathway' and 'OUTPUT-COMPOUNDS' in padmetSpec.dicOfNode[node_pathway].misc]
    if sum(padmet_input_compounds_in_pwys) == 0 or sum(padmet_output_compounds_in_pwys) == 0:
        sys.exit("The padmet " + padmet_file + " does not contain INPUT-COMPOUNDS and OUTPUT-COMPOUNDS in the pathway node, padmet can't produce the pathway visualization without them.")

    edges = []
    edges_label = []
    weights = []

    nodes = {}
    nodes_label = []

    all_pwys = [node for node in padmetSpec.dicOfNode if padmetSpec.dicOfNode[node].type == "pathway"]

    for pwy in all_pwys:
        node_pwy = padmetSpec.dicOfNode[pwy]
        if 'INPUT-COMPOUNDS' in node_pwy.misc and 'OUTPUT-COMPOUNDS' in node_pwy.misc:
            reactants = node_pwy.misc['INPUT-COMPOUNDS'][0].split(',')
            products = node_pwy.misc['OUTPUT-COMPOUNDS'][0].split(',')
            for reactant in reactants:
                if reactant not in nodes:
                    new_cpd_id = len(nodes_label)
                    nodes_label.append(reactant)
                    nodes[reactant] = new_cpd_id
            for product in products:
                if product not in nodes:
                    new_cpd_id = len(nodes_label)
                    nodes_label.append(product)
                    nodes[product] = new_cpd_id
            reactant_product_tuples = [(reactant, product) for reactant in reactants for product in products]
            for reactant, product in reactant_product_tuples:
                edges.append((nodes[reactant], nodes[product]))
                weights.append(1)
                edges_label.append(node_pwy.id)

    return edges, edges_label, weights, nodes, nodes_label


def create_graph(metabolic_network_file, output_file, visualization_level, hide_currency_metabolites):
    """ Using output of parse_compounds_padmet or parse_compounds_sbml create a network picture using igraph.

    Parameters
    ----------
    metabolic_network_file: str
        pathname of the metabolic network file
    output_file: str
        pathname of the output picture of the metabolic network
    visualization_level: str
        level of visualization either compound, reaction or pathway
    hide_currency_metabolites: bool
        hide currency metabolites
    """
    if hide_currency_metabolites:
        hide_metabolites = ["PROTON", "WATER", "OXYGEN-MOLECULE", "NADP", "NADPH", "ATP",
                            "PPI", "CARBON-DIOXIDE", "Pi", "ADP", "CO-A", "UDP", "NAD",
                            "NADH", "AMP", "AMMONIA", "HYDROGEN-PEROXIDE", "Acceptor",
                            "Donor-H2", "3-5-ADP", "GDP", "CARBON-MONOXIDE", "GTP", "FAD"]
    else:
        hide_metabolites = []

    if visualization_level == "compound":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_compounds_padmet(metabolic_network_file, hide_metabolites)
        elif metabolic_network_file.endswith('.sbml'):
            edges, edges_label, weights, nodes, nodes_label = parse_compounds_sbml(metabolic_network_file, hide_metabolites)
        else:
            sys.exit('No correct extension file as input. Input must be a .padmet or a .sbml file.')

    elif visualization_level == "reaction":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_reactions_padmet(metabolic_network_file)
        else:
            sys.exit('No correct extension file as input for pathway level. Input must be a .padmet.')

    elif visualization_level == "pathway":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_pathways_padmet(metabolic_network_file)
        else:
            sys.exit('No correct extension file as input for pathway level. Input must be a .padmet.')

    n_vertices = len(nodes)

    # igraph network implementation is from https://gist.github.com/Vini2/d13b12b37b01b4001cabf38b1f850d8a#file-visualise_graph_demo-ipynb
    # Thanks to Vini2.

    # Create graph
    compounds_graph = igraph.Graph(directed=True)

    # Add vertices
    compounds_graph.add_vertices(n_vertices)

    # Add edges to the graph
    compounds_graph.add_edges(edges)

    # Nodes label
    compounds_graph.vs["label"] = nodes_label

    # Add weights to edges in the graph
    compounds_graph.es['weight'] = weights

    visual_style = {}

    # Define colors used for outdegree visualization
    colours = ['#fecc5c', '#a31a1c']

    # Set bbox and margin
    visual_style["bbox"] = (5000,5000)
    visual_style["margin"] = 17

    # Set vertex colours
    visual_style["vertex_color"] = 'grey'

    # Set vertex size
    visual_style["vertex_size"] = 20

    # Set vertex lable size
    visual_style["vertex_label_size"] = 8

    # Don't curve the edges
    visual_style["edge_curved"] = True

    visual_style["edge_arrow_size"] = 0.5

    # Set the layout
    my_layout = compounds_graph.layout_fruchterman_reingold()
    visual_style["layout"] = my_layout

    # Plot the graph
    igraph.plot(compounds_graph, output_file, **visual_style)


def create_html_graph(metabolic_network_file, output_file, visualization_level):
    """ Using output of parse_compounds_padmet or parse_compounds_sbml create an interactive graph in html.

    Parameters
    ----------
    metabolic_network_file: str
        pathname of the metabolic network file
    output_file: str
        pathname of the output picture of the metabolic network
    visualization_level: str
        level of visualization either compound, reaction or pathway
    """
    try:
        import pyvis
    except ImportError:
        raise ImportError('Requires pyvis, try:\npip install pyvis')

    if visualization_level == "compound":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_compounds_padmet(metabolic_network_file)
        elif metabolic_network_file.endswith('.sbml'):
            edges, edges_label, weights, nodes, nodes_label = parse_compounds_sbml(metabolic_network_file)
        else:
            sys.exit('No correct extension file as input. Input must be a .padmet or a .sbml file.')

    elif visualization_level == "reaction":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_reactions_padmet(metabolic_network_file)
        else:
            sys.exit('No correct extension file as input for pathway level. Input must be a .padmet.')

    elif visualization_level == "pathway":
        if metabolic_network_file.endswith('.padmet'):
            edges, edges_label, weights, nodes, nodes_label = parse_pathways_padmet(metabolic_network_file)
        else:
            sys.exit('No correct extension file as input for pathway level. Input must be a .padmet.')
    
    compounds_graph = pyvis.network.Network(height=1900, width=1900, directed=True)
    compounds_graph.barnes_hut()

    compounds_graph.add_nodes([node for node in nodes.values()], label=[nodes_label[node] for node in nodes.values()])

    for edge in edges:
        compounds_graph.add_edge(edge[0], edge[1])

    compounds_graph.show(output_file)
