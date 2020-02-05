import cobra
import csv
import libsbml
import os
import shutil

from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet.utils.connection.sbmlGenerator import padmet_to_sbml
from padmet.utils.exploration.compare_padmet import compare_padmet
from padmet.utils.exploration.padmet_stats import compute_stats
from padmet.utils.exploration.get_pwy_from_rxn import extract_pwys
from padmet.utils.exploration.flux_analysis import flux_analysis

FABO_RXNS = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN',
            'ENOYL-COA-DELTA-ISOM-RXN', 'OHBUTYRYL-COA-EPIM-RXN', 'KETOACYLCOATHIOL-RXN',
            'OHACYL-COA-DEHYDROG-RXN']

FABO_GENES = ['b2341', 'b0221', 'b1805', 'b1701', 'b2342', 'b3845', 'b3846']

FABO_CPDS = ['D-3-HYDROXYACYL-COA', 'L-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA',
            'CIS-DELTA3-ENOYL-COA', 'Saturated-Fatty-Acyl-CoA', 'ETF-Oxidized',
            'ETF-Reduced', 'WATER', 'PROTON', 'CPD66-39', 'CO-A', 'ATP', 'PPI', 'AMP', 'ACETYL-COA',
            '3-KETOACYL-COA', 'NAD', 'NADH']

def test_compare_padmet():
    fabo_1_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_1_padmetSpec.delNode('ACYLCOASYN-RXN')
    fabo_1_padmetSpec.generateFile('fabo_1.padmet')

    fabo_2_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_2_padmetSpec.delNode('ACYLCOADEHYDROG-RXN')
    fabo_2_padmetSpec.generateFile('fabo_2.padmet')

    compare_padmet('fabo_1.padmet,fabo_2.padmet', 'output', padmetRef = None, verbose = False)

    genes_fabo_1 = []
    genes_fabo_2 = []
    with open('output/genes.csv', 'r') as genes_file:
        csvreader = csv.reader(genes_file, delimiter='\t')
        for row in csvreader:
            if row[1] == 'present':
                genes_fabo_1.append(row[0])
            if row[2] == 'present':
                genes_fabo_2.append(row[0])

    assert set(FABO_GENES).issubset(set(genes_fabo_1))

    assert set(FABO_GENES).issubset(set(genes_fabo_2))

    reactions_fabo_1 = []
    reactions_fabo_2 = []
    with open('output/reactions.csv', 'r') as reactions_file:
        csvreader = csv.reader(reactions_file, delimiter='\t')
        for row in csvreader:
            if row[1] == 'present':
                reactions_fabo_1.append(row[0])
            if row[2] == 'present':
                reactions_fabo_2.append(row[0])

    expected_fabo_1_rxns = [rxn for rxn in FABO_RXNS if rxn != 'ACYLCOASYN-RXN']
    expected_fabo_2_rxns = [rxn for rxn in FABO_RXNS if rxn != 'ACYLCOADEHYDROG-RXN']

    assert set(expected_fabo_1_rxns).issubset(set(reactions_fabo_1))

    assert set(expected_fabo_2_rxns).issubset(set(reactions_fabo_2))

    pathway_fabo_1 = []
    pathway_fabo_2 = []
    with open('output/pathways.csv', 'r') as pathways_file:
        csvreader = csv.reader(pathways_file, delimiter='\t')
        for row in csvreader:
            if row[0] != 'pathway':
                pathway_fabo_1.append(row[0])
                pathway_fabo_2.append(row[0])
            if row[3] != 'fabo_1_rxn_assoc (sep=;)':
                pwy_reactions_fabo_1 = row[3].split(';')
            if row[4] != 'fabo_2_rxn_assoc (sep=;)':
                pwy_reactions_fabo_2 = row[4].split(';')

    assert pathway_fabo_1 == ['FAO-PWY']

    assert pathway_fabo_2 == ['FAO-PWY']

    assert set(expected_fabo_1_rxns).issubset(set(pwy_reactions_fabo_1))

    assert set(expected_fabo_2_rxns).issubset(set(pwy_reactions_fabo_2))


    metabolites_fabo_1 = []
    metabolites_fabo_2 = []
    with open('output/metabolites.csv', 'r') as metabolites_file:
        csvreader = csv.reader(metabolites_file, delimiter='\t')
        for row in csvreader:
            if row[1] != 'fabo_1_rxn_consume' or row[1] != '':
                if row[0] != 'metabolite':
                    metabolites_fabo_1.append(row[0])
            if row[3] != 'fabo_1_rxn_produce' or row[3] != '':
                if row[0] != 'metabolite':
                    metabolites_fabo_1.append(row[0])
            if row[2] != 'fabo_2_rxn_consume' or row[2] != '':
                if row[0] != 'metabolite':
                    metabolites_fabo_2.append(row[0])
            if row[2] != 'fabo_2_rxn_produce' or row[2] != '':
                if row[0] != 'metabolite':
                    metabolites_fabo_2.append(row[0])

    metabolites_fabo_1 = list(set(metabolites_fabo_1))
    metabolites_fabo_2 = list(set(metabolites_fabo_2))

    assert set(FABO_CPDS).issubset(set(metabolites_fabo_1))

    assert set(FABO_CPDS).issubset(set(metabolites_fabo_2))

    os.remove('fabo_1.padmet')
    os.remove('fabo_2.padmet')
    shutil.rmtree('output')


def test_padmet_stats():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_padmetSpec.generateFile('fabo.padmet')
    compute_stats('fabo.padmet')

    # Expected stats: nb pathways, nb reactions, nb reactions with gene, nb genes, nb compounds
    expected_stats = ['1', '10', '7', '8', '25']
    with open('padmet_stats.tsv', 'r') as reactions_file:
        csvreader = csv.reader(reactions_file, delimiter='\t')
        for row in csvreader:
            if row[0] != 'padmet_file':
                fabo_stats = [row[1], row[2], row[3], row[4], row[5]]

    assert fabo_stats == expected_stats

    os.remove('fabo.padmet')
    os.remove('padmet_stats.tsv')
    os.remove('padmet_orthologs_stats.tsv')


def test_get_pwy_from_rxn():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    pathways = extract_pwys(fabo_padmetSpec, set(FABO_RXNS))

    assert pathways['FAO-PWY']['ratio'] == 1.0

    assert pathways['FAO-PWY']['total_rxn'] == set(FABO_RXNS)

    assert pathways['FAO-PWY']['rxn_from_list'] == set(FABO_RXNS)


def test_flux_analysis():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml', obj_fct='KETOACYLCOATHIOL-RXN')
    flux_analysis('fabo.sbml', seeds_file='test_data/seeds.sbml', targets_file='test_data/targets.sbml', all_species=True)
    os.remove('fabo.sbml')
test_flux_analysis()