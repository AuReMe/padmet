# -*- coding: utf-8 -*-
"""
Description:
    Visualize similarity between metabolic networks using MDS.

::

    usage:
        padmet visu_similarity_gsmn --reaction=FILE --group=FILE  --output=FILE
    
    options:
        -h --help     Show help.
        --reaction=FILE    pathname to the reaction file output of compare_padmet or compare_sbml.
        --group=FILE    pathname to the group file containing a column named "species" with the organism ID and a column "group" classifying species in group
        --output=FILE    pathname to the picture output file containing the MDS projection
"""
import docopt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pa
import matplotlib.pyplot as plt

from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.manifold import MDS
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.feature_selection import SelectFromModel


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def visu_similarity_gsmn_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    reaction_file = args["--reaction"]
    group_file = args["--group"]
    output_file = args["--output"]
    visu_similarity_gsmn(reaction_file, group_file, output_file)


def visu_similarity_gsmn(reaction_file, group_file, output_file):
    """
    Create dendrogram, upset figure (if upset argument) and compare reactiosn in species.

    Parameters
    ----------
    reaction_file: str
        path to reaction file from compare_padmet/compare_sbml.
    group_file: str
        path to group file containing group assignation for each metabolic network.
    output_file: str
        path to picture ouput file.
    """
    # Convert the reaction table into a matrix.
    df = pa.read_csv(reaction_file, sep='\t')
    df.set_index('reaction', inplace=True)
    df = df.replace('present', 1)
    df = df.replace(np.nan, 0)

    df = df.transpose()

    # Retrieve the group from the group file.
    df_family = pa.read_csv(group_file, sep='\t')
    df_family.set_index('species', inplace=True)
    df_join = df.join(df_family)

    groups = {index: group for index, group in enumerate(df_join['group'].unique().tolist())}

    # Select only species with a group.
    df_join = df_join[df_join['group'].isin([groups[group] for group in groups])]

    for group in groups:
        df_join = df_join.replace(groups[group], group)

    X = df_join[df.columns.tolist()]
    y = df_join['group']

    # Projection with MDS.
    embedding = MDS(n_components=2)
    X_transformed = embedding.fit_transform(X)
    X = X_transformed

    # Create color map for each group in the figure.
    color_map = plt.cm.get_cmap('hsv', len(groups))

    plt.rcParams['figure.figsize'] = [20, 20]
    plt.rc('font', size=14)

    # Plot each point in a matplotlib figure. 
    for i in np.unique(df_join.group):
        group_color = np.asarray(color_map(i)).reshape(1,-1)
        subset = X[df_join.group == i]
        labels = df_join[df_join.group == i].index
        
        x = [row[0] for row in subset]
        y = [row[1] for row in subset]
        plt.scatter(x,y,c=group_color,label=groups[i])
        for i in range(len(x)):
            plt.annotate(labels[i], (x[i],y[i]))
        plt.legend()

    plt.savefig(output_file)
