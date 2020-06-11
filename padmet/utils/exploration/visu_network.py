# -*- coding: utf-8 -*-
"""
Description:
    Allows to visualize a metabolic network on a compounds perspectives

::

    usage:
        padmet visu_network -i=FILE -o=FILE [--html=FILE]
    
    options:
        -h --help     Show help.
        -i=FILE    pathname to the input file (either PADMet or SBML).
        -o=FILE    pathname to the output file (picture of metabolic network)
        --html=FILE    pathname to the output file (interactive hmtl of metabolic network)
"""
import docopt
import logging

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
    create_compounds_graph(metabolic_network_file, output_file)
    if html_output_file:
        create_html_compounds_graph(metabolic_network_file, html_output_file)


def parse_compounds_padmet(padmet_file):
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
            for compound_out in outs:
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

    return edges, edges_label, weights, nodes, nodes_label


def parse_compounds_sbml(sbml_file):
    sbml_model = read_sbml_model(sbml_file)

    edges = []
    edges_label = []
    weights = []

    nodes = {}
    nodes_label = []

    for reaction in sbml_model.reactions:
        for reactant in reaction.reactants:
            reactant = convert_from_coded_id(reactant.id)[0]
            for product in reaction.products:
                product = convert_from_coded_id(product.id)[0]
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

    return edges, edges_label, weights, nodes, nodes_label


def create_compounds_graph(metabolic_network_file, output_file):

    if metabolic_network_file.endswith('.padmet'):
        edges, edges_label, weights, nodes, nodes_label = parse_compounds_padmet(metabolic_network_file)
    elif metabolic_network_file.endswith('.sbml'):
        edges, edges_label, weights, nodes, nodes_label = parse_compounds_sbml(metabolic_network_file)

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


def create_html_compounds_graph(metabolic_network_file, output_file):
    try:
        import pyvis
    except ImportError:
        raise ImportError('Requires pyvis, try:\npip install pyvis')

    if metabolic_network_file.endswith('.padmet'):
        edges, edges_label, weights, nodes, nodes_label = parse_compounds_padmet(metabolic_network_file)
    elif metabolic_network_file.endswith('.sbml'):
        edges, edges_label, weights, nodes, nodes_label = parse_compounds_sbml(metabolic_network_file)

    
    G = pyvis.network.Network(height=1900, width=1900, directed=True)
    G.barnes_hut()

    G.add_nodes([node for node in nodes.values()], label=[nodes_label[node] for node in nodes.values()])

    for edge in edges:
        G.add_edge(edge[0], edge[1])

    G.show(output_file)


if __name__ == "__main__":
    main()

