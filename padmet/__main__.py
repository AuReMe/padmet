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
    biggAPI_to_padmet    Create Padmet from BIGG API
    check_orthology_input    Check proteomes for AuReMe
    check_orthology_input    Check proteomes for AuReMe
    enhanced_meneco_output    Extract meneco output and add reaction informations
    extract_orthofinder    Extract result from OrthoFinder and creates SBML outputs
    extract_rxn_with_gene_assoc    Extract from a SBML the reactions with gene association
    gbk_to_faa    Extract the proteome from a genbank file
    gene_to_targets    From a list of genes and a PADmet extract the targets linked to each gene
    get_metacyc_ontology    From a PADmet Ref file of MetaCyc create an XML ontology file
    modelSeed_to_padmet    Create a PADmet file from modelSeed reactions and pathways files
    padmet_to_asp    Convert PADmet to ASP
    padmet_to_matrix    Create a a stoichiometry matrix from a PADmet file
    padmet_to_padmet    Merge multiple PADmets into one PADmet file
    padmet_to_tsv    Convert a PADmet file into AskOmics compatible tsv files
    pgdb_to_padmet    Create PADmet from dat files
    sbml_to_curation    Extract reaction(s) from SBML file into an AuReMe compatible form
    sbml_to_padmet    Convert a SBML file into a PADmet file
    sbml_to_sbml    Convert a SBML file into another SBML (level change)
    sbmlGenerator    Create a SBML file using a PADmet file
    wikiGenerator    Create Wiki files using a PADmet file

exploration:
    compare_padmet    Compare n PADmet files
    compare_sbml    Compare n SBML files
    compare_sbml_padmet    Compare a SBML file with a PADmet file
    convert_sbml_db    Use MetaNetX to convert a SBML from one database to another database
    dendrogram_reactions_distance    Create a reaction absence/presence dendrogram using reactions.tsv file (from compare_padmet)
    flux_analysis    Run Flux Balance Analysis on a defined reaction
    get_pwy_from_rxn    Return Pathway using a list of reactions from a PADmet Ref file
    padmet_stats    Compute stats on a PADmet file
    prot2genome    Perfoms the structural annotation of AuCoMe
    report_network    Creates reports for a PADmet file
    visu_path    Allows to visualize a pathway in PADmet network

management:
    manual_curation    Update a PADmet Spec by filling specific forms
    padmet_compart    For a given padmet file, check and update compartment
    padmet_medium    For a given set of compounds, create 2 reactions for each compounds to maintain consistency of the network for flux analysis.
    relation_curation    Manually curate PADmet relation

See 'aucome <command> -h' for more information on a specific command.
"""
import docopt
import sys

from padmet.utils.connection import pgdb_to_padmet
from padmet.utils.connection import biggAPI_to_padmet
from padmet.utils.connection import check_orthology_input
from padmet.utils.connection import enhanced_meneco_output
from padmet.utils.connection import extract_orthofinder
from padmet.utils.connection import extract_rxn_with_gene_assoc
from padmet.utils.connection import gbk_to_faa
from padmet.utils.connection import gene_to_targets
from padmet.utils.connection import get_metacyc_ontology
from padmet.utils.connection import modelSeed_to_padmet
from padmet.utils.connection import padmet_to_asp
from padmet.utils.connection import padmet_to_matrix
from padmet.utils.connection import padmet_to_padmet
from padmet.utils.connection import padmet_to_tsv
from padmet.utils.connection import sbml_to_curation_form
from padmet.utils.connection import sbml_to_padmet
from padmet.utils.connection import sbml_to_sbml
from padmet.utils.connection import sbmlGenerator
from padmet.utils.connection import wikiGenerator

from padmet.utils.exploration import compare_padmet
from padmet.utils.exploration import compare_sbml
from padmet.utils.exploration import compare_sbml_padmet
from padmet.utils.exploration import convert_sbml_db
from padmet.utils.exploration import dendrogram_reactions_distance
from padmet.utils.exploration import flux_analysis
from padmet.utils.exploration import get_pwy_from_rxn
from padmet.utils.exploration import padmet_stats
from padmet.utils.exploration import prot2genome
from padmet.utils.exploration import visu_path
from padmet.utils.exploration import report_network

from padmet.utils.management import manual_curation
from padmet.utils.management import padmet_compart
from padmet.utils.management import padmet_medium
from padmet.utils.management import relation_curation

from padmet.utils import gbr

def main(args=None):
    args = docopt.docopt(__doc__, options_first=True)

    command = args.pop('<command>')
    command_args = args.pop('<args>')

    commands_path = {
        'gbr': gbr,
        'biggAPI_to_padmet': biggAPI_to_padmet,
        'check_orthology_input': check_orthology_input,
        'enhanced_meneco_output': enhanced_meneco_output,
        'extract_orthofinder': extract_orthofinder,
        'extract_rxn_with_gene_assoc': extract_rxn_with_gene_assoc,
        'gbk_to_faa': gbk_to_faa,
        'gene_to_targets': gene_to_targets,
        'get_metacyc_ontology': get_metacyc_ontology,
        'modelSeed_to_padmet': modelSeed_to_padmet,
        'padmet_to_asp': padmet_to_asp,
        'padmet_to_matrix': padmet_to_matrix,
        'padmet_to_padmet': padmet_to_padmet,
        'padmet_to_tsv': padmet_to_tsv,
        'pgdb_to_padmet': pgdb_to_padmet,
        'sbml_to_curation_form': sbml_to_curation_form,
        'sbml_to_padmet': sbml_to_padmet,
        'sbml_to_sbml': sbml_to_sbml,
        'sbmlGenerator': sbmlGenerator,
        'wikiGenerator': wikiGenerator,
        'compare_padmet': compare_padmet,
        'compare_sbml': compare_sbml,
        'compare_sbml_padmet': compare_sbml_padmet,
        'convert_sbml_db': convert_sbml_db,
        'dendrogram_reactions_distance': dendrogram_reactions_distance,
        'flux_analysis': flux_analysis,
        'get_pwy_from_rxn': get_pwy_from_rxn,
        'padmet_stats': padmet_stats,
        'prot2genome': prot2genome,
        'visu_path': visu_path,
        'report_network': report_network,
        'manual_curation': manual_curation,
        'padmet_compart': padmet_compart,
        'padmet_medium': padmet_medium,
        'relation_curation': relation_curation,
    }

    if command:
        if command not in commands_path:
            sys.exit(command + ' not a valid command.')

        if '-h' in command_args:
            # Return help for the command
            getattr(commands_path[command], 'command_help')()
            sys.exit()

        # Add command to command_args to be parse by docopt.
        command_args.insert(0,command)

        # Call the Command-line function of the script
        getattr(commands_path[command], command + '_cli')(command_args)


if __name__ == "__main__":
    main()