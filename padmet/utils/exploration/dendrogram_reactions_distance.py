#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
    Use reactions.tsv file from compare_padmet.py to create a dendrogram using a Jaccard distance.
    
    From the matrix absence/presence of reactions in different species computes a Jaccard distance between these species.
    Apply a hierarchical clustering on these data with a complete linkage. Then create a dendrogram.
    Apply also intervene to create an upset graph on the data.

::

    usage:
        padmet dendrogram_reactions_distance --reactions=FILE --output=FOLDER [--padmetRef=STR] [--pvclust] [--upset=INT] [-v]

    option:
        -h --help    Show help.
        --reactions=FILE    pathname of the file containing reactions in each species of the comparison.
        --output=FOLDER    path to the output folder.
        --pvclust    launch pvclust dendrogram using R
        --padmetRef=STR    path to the padmet Ref file
        -u --upset=INT    number of cluster in the upset graph.
        -v    verbose mode.
"""

import csv
import docopt
import pandas as pa
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import subprocess

sns.set_style("white")
sns.set('poster', rc={'figure.figsize':(150,140), 'lines.linewidth': 10}, font_scale=4)

from collections import defaultdict
from lxml import etree
from padmet.classes import PadmetRef
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage, to_tree
from scipy.spatial.distance import pdist, squareform
from supervenn import supervenn


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def dendrogram_reactions_distance_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    reaction_pathname = args['--reactions']
    upset_cluster = int(args['--upset']) if args['--upset'] else None
    output_pathname = args['--output']
    padmet_ref_file = args['--padmetRef']
    pvclust = args['--pvclust']
    #verbose = args['-v']

    reaction_figure_creation(reaction_pathname, output_pathname, upset_cluster, padmet_ref_file, pvclust)


def pvclust_dendrogram(reactions_dataframe, organisms, output_folder):
    """
    Using a distance matrix, pvclust R package (with rpy2 package) create a dendrogram with bootstrap values.

    Parameters
    ----------
    reactions_dataframe: pandas.DataFrame
        Reactions absence/presence matrix
    organisms: list
        organisms names
    output_folder: str
        path to the output folder

    """
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri

    pvclust = importr("pvclust")
    grdevices = importr('grDevices')
    ape = importr('ape')

    # Make pandas dataframe compatible with R dataframe.
    pandas2ri.activate()

    # Launch pvclust on the data silently and in parallel.
    result = pvclust.pvclust(reactions_dataframe, method_dist="binary", method_hclust="complete", nboot=10000, quiet=True, parallel=True)

    # Create the dendrogram picture.
    grdevices.png(file=output_folder+"/"+"pvclust_reaction_dendrogram.png", width=2048, height=2048, pointsize=24)
    pvclust.plot_pvclust(result)
    grdevices.dev_off()

    # Dendrogram to newick
    hclust_result = result.rx2("hclust")
    phylo_result = ape.as_phylo(hclust_result)
    ape.write_tree(phylo_result, file=output_folder+"/"+"dendrogram.nwk", tree_names=True, digits=2)


def hclust_to_xml(linkage_matrix):
    """
    Using a distance matrix from scipy linkage, create a xml tree corresponding to the hierarchical clustering. Return the root of the tree.

    Parameters
    ----------
    linkage_matrix: ndarray
        linkage matrix

    Returns
    -------
    root:
        root of the xml tree
    """
    _, node_list = to_tree(linkage_matrix, rd=True)
    len_longest_cluster_id = len(str(max([node.id for node in node_list])))
    parent_nodes = {}
    # Create an xml tree from the dendrogram with all the node.
    # Begin by the last cluster (containing all the other cluster).
    # Check child of the cluster to keep the hierarchy.
    for index, node in enumerate(reversed(node_list)):
        if index == 0:
            root = etree.Element('cluster_' + str(node.id).zfill(len_longest_cluster_id))
            if node.get_left():
                parent_nodes[node.get_left().id] = root
            if node.get_right():
                parent_nodes[node.get_right().id] = root
        else:
            subroot = etree.SubElement(parent_nodes[node.id], 'cluster_' + str(node.id).zfill(len_longest_cluster_id))
            if node.get_left():
                parent_nodes[node.get_left().id] = subroot
            if node.get_right():
                parent_nodes[node.get_right().id] = subroot

    return root


def create_intersection_files(root, cluster_leaf_species, reactions_dataframe, output_folder_tree_cluster, metacyc_to_ecs):
    """
    Create intersection files.

    Parameters
    ----------
    root: root
        root of the xml tree
    cluster_leaf_species: dictionary
        for each leaf give the organisms in it
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    output_folder_tree_cluster: str
        path to the output folder
    metacyc_to_ecs: dictionary
        mapping of metayc reaction to EC number

    Returns
    -------
    reactions_clust: dictionary
        reactions in each cluster of the tree
    """
    # Extract reactions using XML Tree.
    cluster_folder_name = {}
    reactions_clust = {}
    for element in root.iter():
        folder_name = element.tag
        intersect_ancestor_reactions = []
        ancestors = []
        for ancestor in element.iterancestors():
            ancestor_reactions = reactions_dataframe[reactions_dataframe[cluster_leaf_species[ancestor.tag]].all(1)==True]
            intersect_ancestor_reactions.extend(ancestor_reactions.index.tolist())
            ancestors.append(ancestor.tag)
        element_reactions = reactions_dataframe[reactions_dataframe[cluster_leaf_species[element.tag]].all(1)==True]
        # Option to select only reactions present in our subgroup and not in other species.
        #element_reactions = element_reactions[element_reactions[list(set(reactions_dataframe.columns.tolist()) - set(cluster_leaf_species[element.tag]))].any(1)==False]
        intersect_element_reactions = element_reactions.index.tolist()
        only_intersect_element = list(set(intersect_element_reactions) - set(intersect_ancestor_reactions))
        for reaction in only_intersect_element:
            if reaction not in cluster_folder_name:
                cluster_folder_name[reaction] = [folder_name]
            else:
                cluster_folder_name[reaction].append(folder_name)
        tmp_reactions_dataframe = reactions_dataframe.loc[only_intersect_element]
        folder_path = '/'.join(reversed(ancestors))
        if not os.path.isdir(output_folder_tree_cluster + folder_path + '/' + element.tag):
            os.mkdir(output_folder_tree_cluster + folder_path + '/' + element.tag)
        if metacyc_to_ecs:
            tmp_reactions_dataframe['EC'] = [','.join(metacyc_to_ecs[reaction]) if reaction in metacyc_to_ecs else np.nan for reaction in tmp_reactions_dataframe.index]
        tmp_reactions_dataframe.to_csv(output_folder_tree_cluster  + folder_path + '/' + folder_name + '/' + folder_name + '.tsv', sep='\t')
        reactions_clust[element.tag] = tmp_reactions_dataframe.index.tolist()

    reactions_dataframe['cluster'] = [','.join(cluster_folder_name[reaction]) if reaction in cluster_folder_name else 'no_cluster' for reaction in reactions_dataframe.index]
    if metacyc_to_ecs:
        reactions_dataframe['EC'] = [','.join(metacyc_to_ecs[reaction]) if reaction in metacyc_to_ecs else np.nan for reaction in reactions_dataframe.index]
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column].replace(True, 1, inplace=True)
        reactions_dataframe[column].replace(False, 0, inplace=True)
    reactions_dataframe.to_csv(output_folder_tree_cluster + 'reaction_cluster.tsv', sep='\t')
    
    return reactions_clust


def create_cluster(reactions_dataframe, absence_presence_matrix, linkage_matrix):
    """
    Cut the dendrogram to create clusters.

    Parameters
    ----------
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    absence_presence_matrix: pandas.DataFrame
        transposition of the reactions dataframe
    linkage_matrix: ndarray
        linkage matrix

    Returns
    -------
    dendrogram_fclusters: dictionary
        {number used to split the linkage matrix: ndarray with the corresponding clusters}
    """
    species_number = len(reactions_dataframe.columns)
    # Extract Dendrogram information using fcluster.
    dendrogram_fclusters = {}
    for i in range(species_number):
        results = fcluster(linkage_matrix, i, criterion='maxclust')
        dendrogram_fclusters[i] = results

    return dendrogram_fclusters


def create_supervenn(absence_presence_matrix, reactions_dataframe, output_folder_upset, dendrogram_fclusters, k, verbose=False):
    """
    Create an supervenn graph.

    Parameters
    ----------
    absence_presence_matrix: pandas.DataFrame
        transposition of the reactions dataframe
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    output_folder_upset: str
        path to output folder
    dendrogram_fclusters: dictionary
        {number used to split the linkage matrix: ndarray with the corresponding clusters}
    k: int
        number of cluster to create
    """
    if k < 2:
        print('supervenn needs at least 2 clusters to work.')
        return
    # Extract species in each cluster.
    results = dendrogram_fclusters[k]

    species = absence_presence_matrix.index.tolist()

    cluster_species = dict(zip(species, results))
    cluster_classes = defaultdict(list)

    for key, value in cluster_species.items():
        cluster_classes[value].append(key)

    # For each group, extract the reactions present in its species to create supervenn sets.
    supervenn_sets = []
    supervenn_labels = []
    for cluster in cluster_classes:
        reactions_temp = []
        for species in cluster_classes[cluster]:
            species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == True]
            reactions_temp.extend(species_reactions_dataframe.index.tolist())
        cluster_reactions[cluster] = set(reactions_temp)

    supervenn(supervenn_sets, supervenn_labels, sets_ordering='minimize gaps')
    plt.savefig(output_folder_upset + '/supervenn.png', bbox_inches='tight')
    plt.clf()

    return


def create_intervene_graph(absence_presence_matrix, reactions_dataframe, temp_data_folder, path_to_intervene, output_folder_upset, dendrogram_fclusters, k, verbose=False):
    """
    Create an upset graph. Deprecated function, no we use supervenn look at create_supervenn function.

    Parameters
    ----------
    absence_presence_matrix: pandas.DataFrame
        transposition of the reactions dataframe
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    temp_data_folder: str
        temporary data folder
    path_to_intervene: str
        path to intervene bin
    output_folder_upset: str
        path to output folder
    dendrogram_fclusters: dictionary
        {number used to split the linkage matrix: ndarray with the corresponding clusters}
    k: int
        number of cluster to create

    """
    if k < 2:
        print('intervene needs at least 2 clusters to work.')
        return
    # Extract species in each cluster.
    results = dendrogram_fclusters[k]

    species = absence_presence_matrix.index.tolist()

    cluster_species = dict(zip(species, results))
    cluster_classes = defaultdict(list)

    for key, value in cluster_species.items():
        cluster_classes[value].append(key)

    # Extract reactions in each cluster.
    cluster_reactions = {}
    for cluster in cluster_classes:
        reactions_temp = []
        for species in cluster_classes[cluster]:
            species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == True]
            reactions_temp.extend(species_reactions_dataframe.index.tolist())
        cluster_reactions[cluster] = set(reactions_temp)

    # Create data for creating upset graph using intervene.
    n = 0
    folder_names = {}
    for cluster in cluster_classes:
        cluster_name = 'upset_cluster_' + str(n)
        df = pa.DataFrame({cluster_name: list(cluster_reactions[cluster])})
        df.to_csv(temp_data_folder+'/'+cluster_name+'.tsv', sep='\t', index=None, header=None)
        folder_names[cluster_name] = '_'.join(cluster_classes[cluster])
        n += 1
    df_cluster_name = pa.DataFrame.from_dict(folder_names, orient='index')
    df_cluster_name.reset_index(inplace=True)
    df_cluster_name.columns = ['species', 'cluster']
    df_cluster_name.to_csv(output_folder_upset+'/cluster_name.tsv', sep='\t', index=None)

    cmd = '{0} upset -i  {1}/*.tsv --type list -o {2} --figtype svg'.format(path_to_intervene, temp_data_folder, output_folder_upset)
    if verbose:
        subprocess.call(cmd, shell=True)
    else:
        FNULL = open(os.devnull, 'w')
        subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    return


def add_dendrogram_node_label(reaction_dendrogram, node_list, reactions_clust, len_longest_cluster_id):
    """
    Using cluster nodes, add label and reactions number on each node of the dendrogram.
    This function comes from this answer on stackoverflow: https://stackoverflow.com/a/43519473 

    Parameters
    ----------
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    node_list: list
        cluster nodes
    reactions_clust: dictionary
        reactions in each cluster of the tree
    len_longest_cluster_id: int
        reactions in each cluster of the tree
    """
    # Add intersection information to dendrogram.
    # Extract coords from dendrogram.
    # Get leave coordinates, which are at y == 0
    Xcoords = [item for sublist in reaction_dendrogram['icoord'] for item in sublist]
    Ycoords = [item for sublist in reaction_dendrogram['dcoord'] for item in sublist]
    leave_coords = [(x,y) for x,y in zip(Xcoords,Ycoords) if y==0]

    # Map leave ID and coords.
    # In the dendogram data structure,
    # leave ids are listed in ascending order according to their x-coordinate
    order = np.argsort([x for x,y in leave_coords])
    id_to_coord = dict(zip(reaction_dendrogram['leaves'], [leave_coords[idx] for idx in order]))

    # Map endpoint of each link to coordinates of parent node.
    # From two childs, compute the parent coords.
    # dendrogram['icoord'] contains 4 x coordinates for U-shape segment.
    # dendrogram['dcoord'] contains 4 y coordinates for U-shape segment.
    # Each 4 couples of x[i]y[i] corresponds to a point in the U-shape segment:
    # y axis
    # |      parent
    # | x1y1___|___x2y2
    # |    |       |
    # |    |       |
    # |    |       |
    # | x0y0       x3y3
    # ------------------ x axis
    # Using x1 and x2 we can compute the center of the segment, which corresponds to x coord of the parent node.
    # Using y1 or y2, we have the y coordinates of the parent node.
    children_to_parent_coords = dict()
    for xcoords, ycoords in zip(reaction_dendrogram['icoord'], reaction_dendrogram['dcoord']):
        xparent_coord = (xcoords[1] + xcoords[2]) / 2
        yparent_coord = ycoords[1]
        parent_coord = (xparent_coord, yparent_coord)
        left_coord = (xcoords[0], ycoords[0])
        right_coord = (xcoords[3], ycoords[3])
        children_to_parent_coords[(left_coord, right_coord)] = parent_coord

    if all((coord[1]==0 for coords in list(children_to_parent_coords.keys()) for coord in coords)) and all((coords[1]==0 for coords in list(children_to_parent_coords.values()))):
        return None

    # Create a range from the latest leaves to the higher node.
    ids_left = range(len(reaction_dendrogram['leaves']), len(node_list))

    # Iterate on all the nodes.
    # Using children (leaves), retrieve the coords of parent (nodes).
    # Until all nodes have coords.
    while len(ids_left) > 0:
        for node_id in ids_left:
            node = node_list[node_id]
            if (node.left.id in id_to_coord) and (node.right.id in id_to_coord):
                left_coord = id_to_coord[node.left.id]
                right_coord = id_to_coord[node.right.id]
                id_to_coord[node_id] = children_to_parent_coords[(left_coord, right_coord)]

        ids_left = [node_id for node_id in range(len(node_list)) if not node_id in id_to_coord]

    # For each node, add the corresponding cluster name and the number of reactions.
    for node_id, (x, y) in id_to_coord.items():
        if not node_list[node_id].is_leaf():
            plt.plot(x, y, 'ro')
            node_label = str(node_id) + ' (' + str(len(reactions_clust['cluster_'+str(node_id).zfill(len_longest_cluster_id)])) + ')'
            plt.annotate(node_label, (x, y), xytext=(0, -8),
                        textcoords='offset points',
                        va='top', ha='center')

    return True

def comparison_cluster(reactions_clust, output_folder_comparison):
    """
    Compare all cluster one against another.

    Parameters
    ----------
    reactions_clust: dictionary
        reactions in each cluster of the tree
    output_folder_comparison: str
        path to output folder
    """
    import itertools
    for cluster_1, cluster_2 in itertools.permutations(reactions_clust, 2):
        test = open(output_folder_comparison + cluster_1 + '_vs_' + cluster_2, 'w')
        test.write(str(set(reactions_clust[cluster_1]) - set(reactions_clust[cluster_2])))
        test.close()


def getNewick(node, newick, parentdist, leaf_names):
    """
    Create a newick file from the root node of the dendrogram.
    This function comes from this answer on stackoverflow: https://stackoverflow.com/a/31878514.

    Parameters
    ----------
    node: scipy.cluster.hierarchy.ClusterNode
        root ClusterNode of the scipy tree
    newick: str
        newick string
    parentdist: str
        root ClusterNode distance from the linkage matrix
    leaf_names: list
        list of organism names
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


def absent_and_specific_reactions(reactions_dataframe, output_folder_tree_cluster, output_folder_specific, output_folder_absent, organisms):
    """
    Compare all cluster one against another.

    Parameters
    ----------
    reactions_dataframe: pandas.DataFrame
        dataframe containing absence/presence of reactions in organism
    output_folder_tree_cluster: str
        path to output tree cluster folder
    output_folder_specific: str
        path to output folder with specific reactions for each species
    output_folder_absent: str
        path to output folder with absent reactions for each species
    organisms: list
        organisms names
    """
    specific_file = output_folder_tree_cluster + 'absent_specific_reactions.tsv'
    specific_output = open(specific_file, 'w')
    specific_writer = csv.writer(specific_output, delimiter='\t')
    specific_writer.writerow(['Organism', 'NB reactions', 'Unique reactions', 'Absent reactions'])
    for species in sorted(organisms):
        reactions_in_species = set(reactions_dataframe[reactions_dataframe[species]==True].index.tolist())
        reactions_absent_in_others = set(reactions_dataframe[reactions_dataframe[list(set(organisms)-{species})].any(1)==False].index.tolist())
        reactions_only_in_species = list(reactions_in_species.intersection(reactions_absent_in_others))
        tmp_reactions_dataframe = reactions_dataframe.loc[reactions_only_in_species]
        tmp_reactions_dataframe.to_csv(output_folder_specific+species+'.tsv', sep='\t')

        reactions_not_in_species = set(reactions_dataframe[reactions_dataframe[species]==False].index.tolist())
        reactions_in_others = set(reactions_dataframe[reactions_dataframe[list(set(organisms)-{species})].all(1)==True].index.tolist())
        reactions_only_not_in_species = list(reactions_not_in_species.intersection(reactions_in_others))
        tmp_reactions_dataframe = reactions_dataframe.loc[reactions_only_not_in_species]
        tmp_reactions_dataframe.to_csv(output_folder_absent+species+'.tsv', sep='\t')
        specific_writer.writerow([species, len(reactions_in_species),
                                        len(list(reactions_in_species.intersection(reactions_absent_in_others))),
                                        len(list(reactions_not_in_species.intersection(reactions_in_others)))])
    specific_output.close()


def create_pvclust_dendrogram(reaction_file, output_folder):
    # Check if output_folder exists, if not create it.
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # Read the reactions file with pandas.
    all_reactions_dataframe = pa.read_csv(reaction_file, sep='\t')
    # Keep column containing absence-presence of reactions.
    # (columns with (sep=;) are column with gene name linked to reactions)
    # (columns with _formula contain the reaction formula)
    columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' not in column]
    columns = [column for column in columns if '_formula' not in column]
    reactions_dataframe = all_reactions_dataframe[columns].copy()

    reactions_dataframe.set_index('reaction', inplace=True)

    # Translate 'present'/(nan) data into a True/False absence-presence matrix.
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column] = [1 if data == 'present' else 0 for data in reactions_dataframe[column]]

    # Extract organisms.
    organisms = reactions_dataframe.index.tolist()

    # Create pvclust dendrogram.
    pvclust_dendrogram(reactions_dataframe, organisms, output_folder)


def reaction_figure_creation(reaction_file, output_folder, upset_cluster=None, padmetRef_file=None, pvclust=None, verbose=False):
    """
    Create dendrogram, upset figure (if upset argument) and compare reactiosn in species.

    Parameters
    ----------
    reaction_file: str
        path to reaction file
    upset_cluster: int
        the number of cluster you want in the intervene figure
    output_folder: str
        path to output folder
    padmet_ref_file: str
        path to padmet ref file
    pvclust: bool
        boolean to launch or not R pvclust dendrogram
    """
    # Check if output_folder exists, if not create it.
    output_folder_tree_cluster = output_folder + '/tree_cluster/'
    output_folder_comparison = output_folder + '/tree_cluster/comparison_cluster/'
    output_folder_specific = output_folder_tree_cluster + 'specific_reactions/'
    output_folder_absent = output_folder_tree_cluster + 'absent_reactions/'
    if upset_cluster:
        output_folder_upset = output_folder + '/upset_graph'
        temp_data_folder = output_folder + '/upset_graph/temp_data/'
        folders = [output_folder, output_folder_tree_cluster, output_folder_comparison, output_folder_specific, output_folder_absent, output_folder_upset, temp_data_folder]
    else:
        folders = [output_folder, output_folder_tree_cluster, output_folder_comparison, output_folder_specific, output_folder_absent]

    for folder in folders:
        if not os.path.isdir(folder):
            os.mkdir(folder)

    path_to_intervene = 'intervene'

    if not os.path.exists(reaction_file):
        raise FileNotFoundError("No reactions.tsv file accessible at " + reaction_file)

    # Read the reactions file with pandas.
    all_reactions_dataframe = pa.read_csv(reaction_file, sep='\t')
    # Keep column containing absence-presence of reactions.
    # (columns with (sep=;) are column with gene name linked to reactions)
    # (columns with _formula contain the reaction formula)
    columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' not in column]
    columns = [column for column in columns if '_formula' not in column]
    reactions_dataframe = all_reactions_dataframe[columns].copy()

    reactions_dataframe.set_index('reaction', inplace=True)

    # Translate 'present'/(nan) data into a True/False absence-presence matrix.
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column] = [True if data == 'present' else False for data in reactions_dataframe[column]]

    # Transpose the matrix to have species as index and reactions as columns.
    absence_presence_matrix = reactions_dataframe.transpose()

    # Compute a distance matrix using the Jaccard distance between species and condense it.
    condensed_distance_matrix_jaccard = pdist(absence_presence_matrix, metric='jaccard')

    # Hierarchical clustering on the condensed distance matrix.
    linkage_matrix = linkage(condensed_distance_matrix_jaccard, method='average', metric='jaccard')

    # Draw a dendrogram of the clustering.
    reaction_dendrogram = dendrogram(linkage_matrix, labels=absence_presence_matrix.index, leaf_font_size=100, leaf_rotation=90)

    # Extract organisms.
    organisms = absence_presence_matrix.index.tolist()

    # Create Newick tree
    tree = to_tree(linkage_matrix,False)
    newick_tree = getNewick(tree, "", tree.dist, organisms)
    newick_path = os.path.join(output_folder,'newick.txt')
    with open(newick_path, 'w') as f:
        f.write(newick_tree)

    # Specific reactions for each species.
    absent_and_specific_reactions(reactions_dataframe, output_folder_tree_cluster, output_folder_specific, output_folder_absent, organisms)

    if pvclust:
        pvclust_reactions_dataframe = all_reactions_dataframe[columns].copy()

        pvclust_reactions_dataframe.set_index('reaction', inplace=True)
        for column in pvclust_reactions_dataframe.columns.tolist():
            pvclust_reactions_dataframe[column] = [1 if data == 'present' else 0 for data in pvclust_reactions_dataframe[column]]
        # Create pvclust dendrogram.
        pvclust_dendrogram(pvclust_reactions_dataframe, organisms, output_folder)

    # Extract all the nodes inside the clustering. 
    _, node_list = to_tree(linkage_matrix, rd=True)

    if padmetRef_file:
        padmet_ref = PadmetRef(padmetRef_file)
        metacyc_to_ecs = {node.id: node.misc['EC-NUMBER'] for node in padmet_ref.dicOfNode.values() if node.type == "reaction" and 'EC-NUMBER' in node.misc}
    else:
        metacyc_to_ecs = {}

    # For each cluster, give the list of organisms in it.
    # Then write it in a file.
    len_longest_cluster_id = len(str(max([node.id for node in node_list])))
    cluster_leaf_species = {}
    for node in node_list:
        node_leafs = node.pre_order(lambda child: organisms[child.id] if child.is_leaf() else None)
        cluster_leaf_species['cluster_'+str(node.id).zfill(len_longest_cluster_id)] = node_leafs

    species_clustered_df = pa.DataFrame(columns=organisms)
    for cluster_leaf in cluster_leaf_species:
        tmp_organism_cluster = [True if organism in cluster_leaf_species[cluster_leaf] else False for organism in species_clustered_df.columns]
        species_clustered_df.loc[cluster_leaf]  = tmp_organism_cluster

    species_clustered_df = species_clustered_df.replace(np.nan, False)
    species_clustered_df.to_csv(output_folder_tree_cluster + 'clustered_species.tsv', sep='\t')

    # Create xml structure from hierarchical clustering.
    root = hclust_to_xml(linkage_matrix)

    # Post order traversal of the tree.
    d = {}
    for element in root.iter():
        d[element.tag] = [child.tag for child in element]

    post_order_clusters = {}
    for node in node_list:
        node_label = 'cluster_'+str(node.id).zfill(len_longest_cluster_id)
        if d[node_label] == []:
            species = cluster_leaf_species[node_label]
            tmp_reactions = reactions_dataframe[reactions_dataframe[species].all(1) == True]
            post_order_clusters[node_label] = tmp_reactions.index.tolist()
        else:
            if set(post_order_clusters[d[node_label][0]]).intersection(set(post_order_clusters[d[node_label][1]])) != set():
                post_order_clusters[node_label] = set(post_order_clusters[d[node_label][0]]).intersection(set(post_order_clusters[d[node_label][1]]))
            else:
                post_order_clusters[node_label] = set(post_order_clusters[d[node_label][0]]).union(set(post_order_clusters[d[node_label][1]]))

    # Use xml structure to create intersection files.
    reactions_clust = create_intersection_files(root, cluster_leaf_species, reactions_dataframe, output_folder_tree_cluster, metacyc_to_ecs)

    comparison_cluster(reactions_clust, output_folder_comparison)

    # Add label contaning cluster name and reaction number to each node.
    check_label = add_dendrogram_node_label(reaction_dendrogram, node_list, reactions_clust, len_longest_cluster_id)

    if not check_label:
        print('Warning: no label for cluster name have been added.')

    # Create dendrogram, bbox option adjsut the figure size.
    plt.savefig(output_folder+'/reaction_dendrogram.png',bbox_inches='tight')
    plt.clf()

    if upset_cluster:
        dendrogram_fclusters = create_cluster(reactions_dataframe, absence_presence_matrix, linkage_matrix)
        create_supervenn(absence_presence_matrix, reactions_dataframe, output_folder_upset, dendrogram_fclusters, k, verbose)
