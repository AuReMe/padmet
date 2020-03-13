import csv
import libsbml
import os
import shutil

from Bio import SeqIO
from padmet.classes.padmetSpec import PadmetSpec
from padmet.utils import sbmlPlugin
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet.utils.connection.sbmlGenerator import padmet_to_sbml
from padmet.utils.connection.sbml_to_padmet import sbml_to_padmetSpec
from padmet.utils.connection.wikiGenerator import wikiGenerator
from padmet.utils.connection.sbml_to_curation_form import sbml_to_curation
from padmet.utils.connection.padmet_to_tsv import padmet_to_tsv
from padmet.utils.connection.gbk_to_faa import gbk_to_faa
from padmet.utils.connection.gene_to_targets import gene_to_targets
from padmet.utils.connection.padmet_to_padmet import padmet_to_padmet
from padmet.utils.connection.padmet_to_matrix import padmet_to_matrix
from padmet.utils.connection.extract_rxn_with_gene_assoc import extract_rxn_with_gene_assoc

FABO_RXNS = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN',
            'ENOYL-COA-DELTA-ISOM-RXN', 'OHBUTYRYL-COA-EPIM-RXN', 'KETOACYLCOATHIOL-RXN',
            'OHACYL-COA-DEHYDROG-RXN']

FABO_PWYS = ['FAO-PWY']

FABO_CPDS = ['D-3-HYDROXYACYL-COA', 'L-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA',
            'CIS-DELTA3-ENOYL-COA', 'Saturated-Fatty-Acyl-CoA', 'ETF-Oxidized',
            'ETF-Reduced', 'WATER', 'PROTON', 'CPD66-39', 'CO-A', 'ATP', 'PPI', 'AMP', 'ACETYL-COA',
            '3-KETOACYL-COA', 'NAD', 'NADH']

FABO_GENES = ['b2341', 'b0221', 'b1805', 'b1701', 'b2342', 'b3845', 'b3846']


def extract_data_padmet(test_padmetSpec):
    total_pwy_id = set()
    total_cpd_id = set()

    all_rxns = [node for node in test_padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node for node in test_padmetSpec.dicOfNode.values() if node.type == "gene"]

    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in test_padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
        pathways_ids = set([rlt.id_out for rlt in test_padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        total_pwy_id.update(pathways_ids)

    all_pwys = [node_id for (node_id, node) in test_padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node_id for (node_id, node) in test_padmetSpec.dicOfNode.items() if node_id in total_cpd_id]
    all_rxns = [node.id for node in test_padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node.id for node in test_padmetSpec.dicOfNode.values() if node.type == "gene"]

    return all_pwys, all_cpds, all_rxns, all_genes


def test_pgdb_to_padmet():

    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    assert test_padmetSpec.dicOfNode['FAO-PWY'].misc['COMMON-NAME'][0] == 'fatty acid &beta;-oxidation I - (generic)'

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))


def test_pgdb_to_padmet_without_genes():
    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb')

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)
    assert all_pwys == FABO_PWYS

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert all_genes == []


def test_pgdb_to_padmet_no_orphan_with_genes():
    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', no_orphan=True, extract_gene=True)

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    # 2.3.1.49-RXN is added as it has no association with genes
    assert '2.3.1.49-RXN' not in all_rxns

    assert set(FABO_CPDS).issubset(set(all_cpds))

    # b24589 is a manually added genes without association to a reaction
    assert 'b24589' not in all_genes


def test_sbmlGenerator():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')
    reader = libsbml.SBMLReader()
    document = reader.readSBML('fabo.sbml')
    model = document.getModel()

    compounds = model.getListOfSpecies()
    reactions = model.getListOfReactions()
    genes = []
    for reactionSBML in reactions:
        notes = sbmlPlugin.parseNotes(reactionSBML)
        if "GENE_ASSOCIATION" in list(notes.keys()):
            # Using sbmlPlugin to recover all genes associated to the reaction
            for gene in sbmlPlugin.parseGeneAssoc(notes["GENE_ASSOCIATION"][0]):
                if gene not in genes:
                    genes.append(gene)
    id_compounds = [sbmlPlugin.convert_from_coded_id(compound.id)[0] for compound in compounds]
    id_reactions = [sbmlPlugin.convert_from_coded_id(reaction.id)[0] for reaction in reactions]

    assert set(FABO_RXNS).issubset(set(id_reactions))

    assert set(FABO_CPDS).issubset(set(id_compounds))

    assert set(FABO_GENES).issubset(set(genes))

    os.remove('fabo.sbml')


def test_sbml_to_padmet():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')
    sbml_to_padmetSpec('fabo.sbml', 'fabo.padmet')
    fabo_padmet = PadmetSpec('fabo.padmet')

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(fabo_padmet)

    assert all_pwys == []

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))

    os.remove('fabo.sbml')
    os.remove('fabo.padmet')


def test_wikiGenerator():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_padmetSpec.generateFile('fabo.padmet')
    wikiGenerator('fabo.padmet', 'output', 'TEST', None, None, None, False)

    os.remove('fabo.padmet')

    test_genes = [gene for gene in os.listdir('output/genes')]

    test_cpds = [metabolite for metabolite in os.listdir('output/metabolites')]

    test_reactions = [reaction for reaction in os.listdir('output/reactions')]

    test_pathways = [pathway for pathway in os.listdir('output/pathways')]

    test_organisms = [organism for organism in os.listdir('output/organisms')]

    test_navigations = [navigation for navigation in os.listdir('output/navigation')]

    assert test_pathways == FABO_PWYS

    assert set(FABO_RXNS).issubset(set(test_reactions))

    assert set(FABO_CPDS).issubset(set(test_cpds))

    assert set(FABO_GENES).issubset(set(test_genes))

    assert test_organisms == ['fabo']

    assert sorted(test_navigations) == sorted(['annotation', 'pathwaytools', 'Category:gene', 'Category:reaction',
                                'MediaWiki:Sidebar', 'Category:pathway', 'Main_Page', 'Category:metabolite',
                                'Category:organism'])

    shutil.rmtree('output')


def test_sbml_to_curation_form():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')
    rxns = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN']
    id_reactions = ['R_'+sbmlPlugin.convert_to_coded_id(reaction) for reaction in rxns]
    sbml_to_curation('fabo.sbml', id_reactions, 'form.txt')
    os.remove('fabo.sbml')

    with open('form.txt', 'r') as form_file:
        form_str = form_file.read()
        for rxn in rxns:
            assert rxn in form_str

    os.remove('form.txt')


def test_gbk_to_fasta():
    gbk_to_faa('test_data/gbk/fatty_acid_beta_oxydation_I_1.gbk', 'fatty_acid_beta_oxydation_I_1.faa', 'locus_tag')

    records = [record for record in SeqIO.parse('fatty_acid_beta_oxydation_I_1.faa', 'fasta')]

    expected_records = [record for record in SeqIO.parse('test_data/gbk/fatty_acid_beta_oxydation_I_1.faa', 'fasta')]

    for index, record in enumerate(records):
        assert record.id == expected_records[index].id
        assert record.seq == expected_records[index].seq

    os.remove('fatty_acid_beta_oxydation_I_1.faa')


def test_gene_to_targets():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    gene_to_targets(fabo_padmetSpec, 'test_data/genes_to_targets.txt', 'targets.txt')

    expected_tagerts = sorted(['PROTON', '3-KETOACYL-COA', 'CIS-DELTA3-ENOYL-COA', 'PPI', 'NADH',
                        'D-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA', 'L-3-HYDROXYACYL-COA',
                        'Saturated-Fatty-Acyl-CoA', 'AMP'])

    with open('targets.txt', 'r') as target_file:
        found_targets = sorted(target_file.read().split('\n'))
    os.remove('targets.txt')
    assert found_targets == expected_tagerts


def test_padmet_to_padmet():
    # Using inpu data, create 2 padmets and delete one reaction in each.
    fabo_1_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_1_padmetSpec.delNode('ACYLCOASYN-RXN')
    fabo_1_padmetSpec.generateFile('fabo_1.padmet')

    _, _, all_rxns, _ = extract_data_padmet(fabo_1_padmetSpec)

    assert not set(FABO_RXNS).issubset(set(all_rxns))

    fabo_2_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    fabo_2_padmetSpec.delNode('ACYLCOADEHYDROG-RXN')
    fabo_2_padmetSpec.generateFile('fabo_2.padmet')

    _, _, all_rxns, _ = extract_data_padmet(fabo_2_padmetSpec)

    assert not set(FABO_RXNS).issubset(set(all_rxns))

    # By merging them we should retrieve the two deleted reactions
    padmet_to_padmet('fabo_1.padmet,fabo_2.padmet', 'fabo.padmet')

    expected_padmet = PadmetSpec('fabo.padmet')

    _, _, all_rxns, _ = extract_data_padmet(expected_padmet)

    assert set(FABO_RXNS).issubset(set(all_rxns))

    os.remove('fabo_1.padmet')
    os.remove('fabo_2.padmet')
    os.remove('fabo.padmet')


def test_padmet_to_padmet_reversibility():
    """
    Test an issue encountered when a reaction is defined in a direction in one padmet.
    And in another the same reaction is defined as reversible.
    In old padmet this leads to the reaction having all of its reactants/products as reactants and also as products.
    """
    # Read padmet file
    padmet_to_padmet('test_data/padmet', 'fabo.padmet')

    expected_padmet = PadmetSpec('fabo.padmet')

    reactants = [rlt.id_out for rlt in expected_padmet.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'] if rlt.type in ['consumes']]
    products = [rlt.id_out for rlt in expected_padmet.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'] if rlt.type in ['produces']]

    assert sorted(reactants) == ['TRANS-D2-ENOYL-COA', 'WATER'] or sorted(products) == ['TRANS-D2-ENOYL-COA', 'WATER']
    assert sorted(products) == ['L-3-HYDROXYACYL-COA'] or sorted(reactants) == ['L-3-HYDROXYACYL-COA']

    os.remove('fabo.padmet')


def test_padmet_to_matrix():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb')
    padmet_to_matrix(fabo_padmetSpec, 'matrix.tsv')

    expected_matrix = []
    with open('test_data/stoechiometry_matrix.tsv', 'r') as expected_output:
        expected_matrix_reader = csv.reader(expected_output, delimiter='\t')
        expected_matrix = [row for row in expected_matrix_reader]

    found_matrix = []
    with open('matrix.tsv', 'r') as found_output:
        found_matrix_reader = csv.reader(found_output, delimiter='\t')
        found_matrix = [row for row in found_matrix_reader]

    assert found_matrix == expected_matrix
    os.remove('matrix.tsv')


def test_extract_rxn_with_gene_assoc():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')

    # Extract reactions with only genes association so 2.3.1.49-RXN should not be here.
    extract_rxn_with_gene_assoc('fabo.sbml', 'fabo_rxn_with_genes.sbml', verbose=False)

    reader = libsbml.SBMLReader()
    document = reader.readSBML('fabo_rxn_with_genes.sbml')
    model = document.getModel()
    reactions = model.getListOfReactions()
    id_reactions = [sbmlPlugin.convert_from_coded_id(reaction.id)[0] for reaction in reactions]

    os.remove('fabo.sbml')
    os.remove('fabo_rxn_with_genes.sbml')

    assert '2.3.1.49-RXN' not in id_reactions
