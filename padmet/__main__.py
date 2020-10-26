"""
usage:
    padmet <command> [<args>...]

options:
    -h --help     Show help.
    -v     Verbose.

The subcommands are:

utils:
    gbr    Extract 'AND' and 'OR' from gene string

connection:
    biggAPI_to_padmet    Create PADMet file from BIGG API
    check_orthology_input    Check proteomes for AuReMe
    check_orthology_input    Check proteomes for AuReMe
    enhanced_meneco_output    Extract meneco output and add reaction informations
    extract_orthofinder    Extract result from OrthoFinder and creates SBML outputs
    extract_rxn_with_gene_assoc    Extract from a SBML the reactions with gene association
    gbk_to_faa    Extract the proteome from a genbank file
    gene_to_targets    From a list of genes and a PADMet extract the targets linked to each gene
    get_metacyc_ontology    From a PADMet Ref file of MetaCyc create an XML ontology file
    modelSeed_to_padmet    Create a PADMet file from modelSeed reactions and pathways files
    padmet_to_asp    Convert PADMet to ASP
    padmet_to_matrix    Create a a stoichiometry matrix from a PADMet file
    padmet_to_padmet    Merge multiple PADMets into one PADMet file
    padmet_to_tsv    Convert a PADMet file into AskOmics compatible tsv files
    pgdb_to_padmet    Create PADMet from dat files
    sbml_to_curation    Extract reaction(s) from SBML file into an AuReMe compatible form
    sbml_to_padmet    Convert a SBML file into a PADMet file
    sbml_to_sbml    Convert a SBML file into another SBML (level change)
    sbmlGenerator    Create a SBML file using a PADMet file
    wikiGenerator    Create Wiki files using a PADMet file

exploration:
    compare_padmet    Compare n PADMet files
    compare_sbml    Compare n SBML files
    compare_sbml_padmet    Compare a SBML file with a PADMet file
    convert_sbml_db    Use MetaNetX to convert a SBML from one database to another database
    dendrogram_reactions_distance    Create a reaction absence/presence dendrogram using reactions.tsv file (from compare_padmet)
    flux_analysis    Run Flux Balance Analysis on a defined reaction
    get_pwy_from_rxn    Return Pathway using a list of reactions from a PADMet Ref file
    padmet_stats    Compute stats on a PADMet file
    prot2genome    Perfoms the structural annotation of AuCoMe
    report_network    Creates reports for a PADMet file
    visu_network    Allows to visualize a PADMet or SBML network
    visu_path    Allows to visualize a pathway in PADMet network
    visu_similarity_gsmn    Visualize similarity between metabolic networks using MDS

management:
    manual_curation    Update a PADMet Spec by filling specific forms
    padmet_compart    For a given PADMet file, check and update compartment
    padmet_medium    For a given set of compounds, create 2 reactions for each compounds to maintain consistency of the network for flux analysis.
    relation_curation    Manually curate PADMet relation

See 'padmet <command> -h' for more information on a specific command.
"""
import docopt
import importlib
import sys


def main(args=None):
    args = docopt.docopt(__doc__, options_first=True)

    command = args.pop('<command>')
    command_args = args.pop('<args>')

    commands_path = {
        'utils': ['gbr'] ,
        'utils.connection': [
            'biggAPI_to_padmet', 'check_orthology_input',
            'enhanced_meneco_output', 'extract_orthofinder',
            'extract_rxn_with_gene_assoc', 'gbk_to_faa',
            'gene_to_targets', 'get_metacyc_ontology',
            'modelSeed_to_padmet', 'padmet_to_asp',
            'padmet_to_matrix', 'padmet_to_padmet',
            'padmet_to_tsv', 'pgdb_to_padmet',
            'sbml_to_curation_form', 'sbml_to_padmet',
            'sbml_to_sbml', 'sbmlGenerator',
            'wikiGenerator',
            ],
        'utils.exploration': [
            'compare_padmet', 'compare_sbml',
            'compare_sbml_padmet', 'convert_sbml_db',
            'dendrogram_reactions_distance', 'flux_analysis',
            'get_pwy_from_rxn', 'padmet_stats',
            'prot2genome', 'visu_path',
            'visu_network', 'visu_similarity_gsmn',
            'report_network',
            ],
        'utils.management': [
            'manual_curation', 'padmet_compart',
            'padmet_medium', 'relation_curation',
        ]
    }

    if command:

        # Find the import path of the module
        command_import_path = None
        for import_path in commands_path:
            if command in commands_path[import_path]:
                command_import_path = import_path

        if not command_import_path:
            sys.exit(command + ' not a valid command.')

        # Import the corresponding module
        command_import = importlib.import_module('.'+command, 'padmet.'+command_import_path)
        if '-h' in command_args:
            # Return help for the command
            getattr(command_import, 'command_help')()
            sys.exit()

        # Add command to command_args to be parse by docopt.
        command_args.insert(0,command)

        # Call the Command-line function of the script
        getattr(command_import, command + '_cli')(command_args)


if __name__ == "__main__":
    main()