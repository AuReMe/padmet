# Changelog

# Padmet 5.0.3 (2024-10-30)

## Fix

* issue in GitHub CI.

# Padmet 5.0.2 (2024-10-30)

## Add

* Add getter for pathway with reactions.
* Add metexploreviz_export script.
* Add script to create input for GNN from padmet or sbml.
* Add pathway_production in padmet exploration.
* Implement adding uniprot ids refs from protein-seq-ids-reduced-70.dat file in pgdb_to_padmet.

## Fix

* Bug in `--enhance` functionality from pgdb_to_padmet, reaction IDs were put in dicOfNode instead of reactant IDs (issue #12). Add a test to check this.
* Check if output file is writable for sbmlGenerator.
* Fix issue with pvclust dendrogram using old reactions.tsv.
* Fix an issue with multiprocessing in prot2genomes.
* Fix numerous issues in visualisation scripts.
* Fix doc building.

## Modification

* Replace `setup.py` and `setup.cfg` by `pyproject.toml`.
* Extract inchi from PGDB.
* Update enhanced_meneco_output.py: add json + reactions options argument to consider json meneco output format and be able to extract other set of reactions than union (intersection, minimal, essential).

# Padmet 5.0.1 (2020-10-26)

## Fix

* An issue in setup.py with PyPI and GitHub Actions.

# Padmet 5.0.0 (2020-10-26)

Padmet has now an entrypoint `padmet` to use all the padmet.utils scripts with the corresponding subcommand. More information in the readme.

## Add

* Continuous Integration by using GitHub Actions.
* New visualization scripts: visu_network.py, visu_path.py and visu_similarity_gsmn.py.
* Script to create ontology from MetaCyc.
* Extraction of protein sequence find by exonerate in prot2genome.
* Option to pgdb_to_padmet to keep self-producing reactions.
* The possibility to create newick file from pvclust dendrogram.
* An option to keep temporary files in prot2genomes.
* Misc in padmet Node about spontaneous reactions to keep them in pgdb_to_padmet.
* Multiprocessing in compare_padmet.
* An option to create a file for from_pgdb_to_padmet in pgdb_to_padmet.
* Input and output compounds to pathway in pgdb_to_padmet.
* Getter for Relation.
* Output folder for padmet_stats.

## Fix

* Dependencies for padmet are now flexible (>=) instead of being fixed (==).
* Consistency in padmet_to_padmet ';' to ','.
* Issue in padmet_to_tsv when there is no padmetRef file.
* Error when trying to create node label in dendrogram_reactions_distance.py.
* Issue with pgdb_to_padmet (complex with compounds/classes/RNAs).
* Issue in pgdb_to_padmet not taking into account attributes with " - ".
* Issue in padmet_to_padmet when merging the same reaction but with different directions.
* Issue with biopython >= 1.78 in gbk_to_faa.
* Numerous typos.

## Modification

* Refactor padmet creation by implementing the instantiate_padmet function and the `padmet/classes/instantiation.py` script.
* Completely rewrite compare_sbml to be able to compare multiple sbml files.
* Better error message for menetools in flux_balance_analysis.
* Update readme and docs.
* Update cobra dependency to 0.17.1 and fix issue with old function of cobra.
* Replace intervene and upset graph by supervenn in dendrogram_reactions_distance.py.
* Rename `.csv` extension into `.tsv` because we only use tab delimiter.
* Regex for SBML ID.
* Replace `present'  in compare_padmet and comapre_sbml by 1 (or 0).

## Optimization

* Use a region - 10000 and + 10000 around the match query of tBlastn to use exonerate on it  instead of using the all contig. 
* Try to reduce the number of Orthologues files analysed in extract_orthofinder.py.

