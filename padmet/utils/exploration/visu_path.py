# -*- coding: utf-8 -*-
"""
Description:
    #TODO: not stable version
    Allows to visualize a pathway in padmet network.
Color code: 
reactions associated to the pathway, present in the network: lightGreen
reactions associated to the pathway, not present in the network: red
compounds: skyblue

::

    usage:
        padmet visu_path --padmetSpec=FILE --padmetRef=FILE --pathway=ID
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE    pathname to the PADMet file of the network.
        --padmetRef=FILE    pathname to the PADMet file of the db of reference.
        --pathway=ID    pathway id to visualize.

"""
import docopt
import matplotlib.pylab as plt
import networkx as nx

from padmet.classes import PadmetRef, PadmetSpec
from networkx.drawing.nx_agraph import graphviz_layout


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def visu_path_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    
    padmet_ref = PadmetRef(args["--padmetRef"])
    padmet = PadmetSpec(args["--padmetSpec"])
    pathway_id = args["--pathway"]
    
    #get all reactions in pathway
    try:
        all_reactions = [rlt.id_in for rlt in padmet_ref.dicOfRelationOut.get(pathway_id,None)
        if rlt.type == "is_in_pathway"]
    except TypeError:
        print("%s not in padmetRef" %pathway_id)
        exit()
    reactions_in_network = []
    for reaction_id in all_reactions:
        try:
            padmet.dicOfNode[reaction_id]
            reactions_in_network.append(reaction_id)
        except KeyError:
            pass
    
    
    DG=nx.DiGraph()
    
    custom_node_color = {}
    for reaction_id in all_reactions:
    
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
            custom_node_color[reac] = "skyblue"
            DG.add_edge(reac, reaction_id)
        for prod in products:
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
                     node_list=list(custom_node_color.keys()),
                     node_color=list(custom_node_color.values()))
    plt.axis('off')
    
    #save_plot(plt, 'pathway_' + pathway_id)
    plt.show()
    
def save_plot(plot, filepath):
    """Saves plot in multiple formats

    :param arg1: plot object
    :param arg2: filename with its path
    :type arg1: <plot object>
    :type arg2: <str>

    """

    plot.savefig(filepath + ".png",
                 dpi=144, format='png')
    plot.savefig(filepath + ".svg",
                 dpi=144, format='svg')
    #plot.savefig(filepath + ".pdf",
    #              dpi=144, format='pdf')


if __name__ == "__main__":
    main()

