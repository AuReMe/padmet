# -*- coding: utf-8 -*-
"""
Description:
    Visualize similarity between metabolic networks using MDS.

::

    usage:
        padmet visu_similarity_gsmn --reaction=FILE --output=FILE [--group=FILE]
    
    options:
        -h --help     Show help.
        --reaction=FILE    pathname to the reaction file output of compare_padmet or compare_sbml.
        --output=FILE    pathname to the picture output file containing the MDS projection
        --group=FILE    pathname to the group file containing a column named "species" with the organism ID and a column "group" classifying species in group (you can also use a "color" column to associate group to specific color)
"""
import docopt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pa

try:
    from sklearn.manifold import MDS
except ImportError:
    raise ImportError('Requires sklearn, try:\npip install scikit-learn')


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def visu_similarity_gsmn_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    reaction_file = args["--reaction"]
    output_file = args["--output"]
    group_file = args["--group"]
    visu_similarity_gsmn(reaction_file, output_file, group_file)


def visu_similarity_gsmn(reaction_file, output_file, group_file=None):
    """
    Create dendrogram, upset figure (if upset argument) and compare reactiosn in species.

    Parameters
    ----------
    reaction_file: str
        path to reaction file from compare_padmet/compare_sbml.
    output_file: str
        path to picture ouput file.
    group_file: str
        path to group file containing group assignation for each metabolic network.
    """
    # Convert the reaction table into a matrix.
    df = pa.read_csv(reaction_file, sep='\t')
    df.set_index('reaction', inplace=True)
    df = df[[column for column in df.columns if "(" not in column and "_formula" not in column]]
    df = df.replace('present', 1)
    df = df.replace(np.nan, 0)

    df = df.transpose()

    if group_file:
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

        if "color" in df_family.columns:
            # Extract color from group file.
            color_map = dict(zip(df_family.group, df_family.color))
        else:
            # Create color map for each group in the figure.
            color_map = plt.cm.get_cmap('hsv', len(groups))
    else:
        X = df

    # Projection with MDS.
    embedding = MDS(n_components=2)
    X_transformed = embedding.fit_transform(X)
    X = X_transformed

    plt.rcParams['figure.figsize'] = [20, 20]
    plt.rc('font', size=14)

    # Plot each point in a matplotlib figure. 
    if group_file:
        # Use group as provided by the user.
        for i in groups:
            if isinstance(color_map, dict):
                group_color = color_map[groups[i]]
            else:
                group_color = np.asarray(color_map(i)).reshape(1,-1)
            subset = X[df_join.group == i]
            labels = df_join[df_join.group == i].index

            x = [row[0] for row in subset]
            y = [row[1] for row in subset]
            plt.scatter(x, y, c=group_color, label=groups[i])
            for i in range(len(x)):
                plt.annotate(labels[i], (x[i],y[i]))
            plt.legend()
    else:
        for genome in df.index:
            subset = X[df.index == genome]

            x = [row[0] for row in subset]
            y = [row[1] for row in subset]
            plt.scatter(x, y, c="black", label=genome)
            for i in range(len(x)):
                plt.annotate(genome, (x[i],y[i]))

    plt.savefig(output_file)
