#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create a padmet stats file containing the number of pathways, reactions,
    genes and compounds inside the padmet.

    The input is a padmet file or a folder containing multiple padmets.
    
    Create a tsv file named padmet_stats.tsv where the script have been
    launched.

::

usage:
    padmet padmet_stats --padmet=FILE --output=FOLDER

option:
    -h --help    Show help.
    -p --padmet=FILE    padmet file or folder containing padmet(s).
    -o --output=FOLDER    path to output folder.
"""

import csv
import docopt
import os

from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def padmet_stats_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    padmet_file_folder = args['--padmet']
    output_folder = args['--output']

    compute_stats(padmet_file_folder, output_folder)


def compute_stats(padmet_file_folder, output_folder):
    """
    Count reactions/pathways/compounds/genes in padmet(s).

    Parameters
    ----------
    padmet_file_folder: str
        path to the padmet file/folder to analyze
    output_folder: str
        path to the output folder
    """
    if os.path.isdir(padmet_file_folder):
        padmet_type = "dir"
    elif os.path.isfile(padmet_file_folder):
        padmet_type = "file"
    else:
        raise TypeError("%s is not a dir or a file or is not accessible." %(padmet_file_folder))

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_file = open(output_folder + '/padmet_stats.tsv', 'w')
    output_writer = csv.writer(output_file, delimiter='\t')
    output_writer.writerow(['padmet_file', 'pathways_wtih_reactions', 'reactions', 'reactions_with_gene_association', 'genes', 'compounds'])


    orthologs_species = []
    orthologs_stats = {}

    if padmet_type == "dir":
        padmet_names = [padmet_file.replace('.padmet', '').upper() for padmet_file in os.listdir(padmet_file_folder)]
        for padmet_file in os.listdir(padmet_file_folder):
            padmet_path = padmet_file_folder + '/' + padmet_file
            stats = padmet_stat(padmet_path)
            output_writer.writerow(stats)
            ortholog_species_counts = orthology_result(padmet_path, padmet_names)
            orthologs_stats[padmet_file.replace('.padmet', '').upper()] = ortholog_species_counts
            for species in ortholog_species_counts:
                if species not in orthologs_species:
                    orthologs_species.append(species)

    if padmet_type == "file":
        stats = padmet_stat(padmet_file_folder)
        output_writer.writerow(stats)
        padmet_name = padmet_file_folder.replace('.padmet', '').upper()
        padmet_names = [padmet_name]
        ortholog_species_counts = orthology_result(padmet_file_folder, padmet_names)
        orthologs_stats[padmet_file_folder.replace('.padmet', '').upper()] = ortholog_species_counts
        for species in ortholog_species_counts:
            if species not in orthologs_species:
                orthologs_species.append(species)

    orthologs_stats_file = open(output_folder + '/padmet_orthologs_stats.tsv', 'w')
    orthologs_stats_writer = csv.writer(orthologs_stats_file, delimiter='\t')
    orthologs_stats_writer.writerow(['', *orthologs_species])
    for input_species in orthologs_stats:
        orthologs_stats_writer.writerow([input_species, *[orthologs_stats[input_species][species] for species in orthologs_species]])

    output_file.close()
    orthologs_stats_file.close()


def padmet_stat(padmet_file):
    """
    Count reactions/pathways/compounds/genes in a padmet file.

    Parameters
    ----------
    padmet_file: str
        path to a padmet file

    Returns
    -------
    list:
        [path to padmet, number of pathways, number of reactions, number of genes, number of compounds]
    """
    padmetSpec = PadmetSpec(padmet_file)
    padmet_name = os.path.basename(padmet_file).replace('.padmet', '')
    total_pwy_id = set()
    total_cpd_id = set()

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node for node in padmetSpec.dicOfNode.values() if node.type == "gene"]

    nb_rxn_with_ga = 0
    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
		# Get all pathways having at least a reaction. Remove superpathways containing only pathways.
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        if any([rlt for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_linked_to"]):
            nb_rxn_with_ga += 1
        total_pwy_id.update(pathways_ids)

    all_pwys = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]

    return [padmet_name, len(all_pwys), len(all_rxns), nb_rxn_with_ga, len(all_genes), len(all_cpds)]


def orthology_result(padmet_file, padmet_names):
    """
    Count reactions/pathways/compounds/genes in a padmet file.

    Parameters
    ----------
    padmet_file: str
        path to a padmet file
    padmet_names: list
        all the padmet filenames

    Returns
    -------
    dictionary:
        Number of reactions given by the other species
    """
    ortholog_species_counts = {}

    padmetSpec = PadmetSpec(padmet_file)

    ortholog_reactions = set()
    for node in padmetSpec.dicOfNode.values():
        if node.type == 'suppData':
            reaction_id = node.id.split('_SuppData_OUTPUT_ORTHOFINDER_')[0]
            ortholog_reactions.add(reaction_id)

            ortholog_species = node.id.split('FROM_')[1]
            if ortholog_species in ortholog_species_counts:
                ortholog_species_counts[ortholog_species] += 1
            else:
                ortholog_species_counts[ortholog_species] = 1

    for species in padmet_names:
        if species not in ortholog_species_counts:
            ortholog_species_counts[species] = 0

    ortholog_species_counts['Orthology'] = len(ortholog_reactions)

    return ortholog_species_counts

