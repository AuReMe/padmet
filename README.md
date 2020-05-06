[![PyPI version](https://img.shields.io/pypi/v/padmet.svg)](https://pypi.org/project/padmet/) [![GitHub license](https://img.shields.io/github/license/AuReMe/padmet.svg)](https://github.com/AuReMe/padmet/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/padmet/badge/?version=latest)](https://padmet.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/AuReMe/padmet.svg?branch=master)](https://travis-ci.org/AuReMe/padmet)

# PADMet

PADMet is a python library to manage metabolic networks.

## Table of contents
- [PADMet](#padmet)
  - [Table of contents](#table-of-contents)
  - [Description](#description)
  - [Installation](#installation)
    - [Installation with pip](#installation-with-pip)
  - [Use](#use)
    - [Entrypoint](#entrypoint)
    - [Python import](#python-import)
    - [Corresponding functions](#corresponding-functions)

## Description

The PADMet package allows conciliating genomics and metabolic network information used to produce a genome-scale constraint-based metabolic model within a database that traces all the reconstruction process steps. It allows representing the metabolic model in the form of a Wiki containing all the used/traced information. Other standard outputs are made available with the package. 

The main concept underlying PADMet-Package is to provide solutions that ensure the consistency, the internal standardization and the reconciliation of the information used within any workflow that combines several tools involving metabolic networks reconstruction or analysis. The PADMet package is at the core of the AuReMe workflow, dedicated to the primary reconstruction of genome-scale metabolic networks from raw data. It allows the study of organisms for which few experimental data are available. Its main feature is to undergo the reconstruction of the metabolic network by combining several heterogeneous knowledge and data sources, including the information reported by several scaffold metabolic networks for cousin species.

## Installation

### Installation with pip

```
pip install padmet
```

## Use

PADMet can be used either by using its entrypoint or by making an import in python.

### Entrypoint

```
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
    visu_path    Allows to visualize a pathway in PADMet network

management:
    manual_curation    Update a PADMet Spec by filling specific forms
    padmet_compart    For a given PADMet file, check and update compartment
    padmet_medium    For a given set of compounds, create 2 reactions for each compounds to maintain consistency of the network for flux analysis.
    relation_curation    Manually curate PADMet relation

See 'padmet <command> -h' for more information on a specific command.
```

### Python import

```python
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet

padmet_instance = from_pgdb_to_padmet(pgdb_folder, extract_gene=True)
```

### Corresponding functions

| script                        | function to import                                                                                | command-line                       |
|-------------------------------|---------------------------------------------------------------------------------------------------|------------------------------------|
| biggAPI_to_padmet             | ```from padmet.utils.connection.biggAPI_to_padmet import biggAPI_to_padmet```                     | padmet biggAPI_to_padmet           |
| check_orthology_input         | ```from padmet.utils.connection.check_orthology_input import check_orthology_input```             | padmet check_orthology_input       |
| enhanced_meneco_output        | ```from padmet.utils.connection.enhanced_meneco_output import check_orthology_input```            | padmet enhanced_meneco_output      |
| extract_orthofinder           | ```from padmet.utils.connection.extract_orthofinder import orthogroups_to_sbml```                 | padmet extract_orthofinder         |
| extract_orthofinder           | ```from padmet.utils.connection.extract_orthofinder import orthologue_to_sbml```                  | padmet extract_orthofinder         |
| extract_rxn_with_gene_assoc   | ```from padmet.utils.connection.extract_rxn_with_gene_assoc import extract_rxn_with_gene_assoc``` | padmet extract_rxn_with_gene_assoc |
| gbk_to_faa                    | ```from padmet.utils.connection.gbk_to_faa import gbk_to_faa```                                   | padmet gbk_to_faa                  |
| gene_to_targets               | ```from padmet.utils.connection.gene_to_targets import gene_to_targets```                         | padmet gene_to_targets             |
| get_metacyc_ontology          | ```from padmet.utils.connection.get_metacyc_ontology import metacyc_to_ontology```                | padmet get_metacyc_ontology        |
| modelSeed_to_padmet           | ```from padmet.utils.connection.modelSeed_to_padmet import modelSeed_to_padmet```                 | padmet modelSeed_to_padmet         |
| padmet_to_asp                 | ```from padmet.utils.connection.padmet_to_asp import padmet_to_asp```                             | padmet padmet_to_asp               |
| padmet_to_matrix              | ```from padmet.utils.connection.padmet_to_matrix import padmet_to_matrix```                       | padmet padmet_to_matrix            |
| padmet_to_padmet              | ```from padmet.utils.connection.padmet_to_padmet import padmet_to_padmet```                       | padmet padmet_to_padmet            |
| padmet_to_tsv                 | ```from padmet.utils.connection.padmet_to_tsv import padmet_to_tsv```                             | padmet padmet_to_tsv               |
| pgdb_to_padmet                | ```from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet```                      | padmet pgdb_to_padmet              |
| sbmlGenerator                 | ```from padmet.utils.connection.sbmlGenerator import padmet_to_sbml```                            | padmet sbmlGenerator               |
| sbml_to_curation_form         | ```from padmet.utils.connection.sbml_to_curation_form import sbml_to_curation```                  | padmet sbml_to_curation_form       |
| sbml_to_padmet                | ```from padmet.utils.connection.sbml_to_padmet import sbml_to_padmetSpec```                       | padmet sbml_to_padmet              |
| sbml_to_sbml                  | ```from padmet.utils.connection.sbml_to_sbml import from_sbml_to_sbml```                          | padmet sbml_to_sbml                |
| wikiGenerator                 | ```from padmet.utils.connection.wikiGenerator import wikiGenerator```                             | padmet wikiGenerator               |
| compare_padmet                | ```from padmet.utils.exploration.compare_padmet import compare_padmet```                          | padmet compare_padmet              |
| compare_sbml                  | ```from padmet.utils.exploration.compare_sbml import compare_sbml```                              | padmet compare_sbml                |
| compare_sbml_padmet           | ```from padmet.utils.exploration.compare_sbml_padmet import compare_sbml_padmet```                | padmet compare_sbml_padmet         |
| convert_sbml_db               | ```from padmet.utils.exploration.convert_sbml_db import map_sbml```                               | padmet convert_sbml_db             |
| dendrogram_reactions_distance | ```from padmet.utils.exploration.dendrogram_reactions_distance import reaction_figure_creation``` |padmet dendrogram_reactions_distance|
| flux_analysis                 | ```from padmet.utils.exploration.flux_analysis import flux_analysis```                            | padmet flux_analysis               |
| get_pwy_from_rxn              | ```from padmet.utils.exploration.get_pwy_from_rxn import get_pwy_from_rxn```                      | padmet get_pwy_from_rxn            |
| padmet_stats                  | ```from padmet.utils.exploration.padmet_stats import compute_stats```                             | padmet padmet_stats                |
| prot2genome                   | ```from padmet.utils.exploration.prot2genome import fromAucome```                                 | padmet prot2genome                 |
| report_network                | ```padmetSpec.network_report(output_dir, padmetRef_file, verbose)```                              | padmet report_network              |
| visu_path                     | ```from padmet.utils.exploration.visu_path import visu_path```                                    | padmet visu_path                   |
| manual_curation               |                                                                                                   | padmet manual_curation             |
| padmet_compart                |                                                                                                   | padmet padmet_compart              |
| padmet_medium                 | ```from padmet.utils.exploration.padmet_medium import manage_medium```                            | padmet padmet_medium               |
| relation_curation             | ```from padmet.utils.exploration.relation_curation import get_relations```                        | padmet relation_curation           |
