import csv
import libsbml
import json
import os
import shutil
import subprocess

from Bio import SeqIO
from padmet.classes.padmetSpec import PadmetSpec
from padmet.utils import sbmlPlugin
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet.utils.connection.sbmlGenerator import padmet_to_sbml
from padmet.utils.connection.sbml_to_sbml import from_sbml_to_sbml
from padmet.utils.connection.sbml_to_padmet import sbml_to_padmetSpec
from padmet.utils.connection.wikiGenerator import wikiGenerator
from padmet.utils.connection.sbml_to_curation_form import sbml_to_curation
from padmet.utils.connection.padmet_to_tsv import padmet_to_tsv
from padmet.utils.connection.gbk_to_faa import gbk_to_faa
from padmet.utils.connection.gene_to_targets import gene_to_targets
from padmet.utils.connection.padmet_to_padmet import padmet_to_padmet
from padmet.utils.connection.padmet_to_matrix import padmet_to_matrix
from padmet.utils.connection.extract_rxn_with_gene_assoc import extract_rxn_with_gene_assoc
from padmet.utils.connection.metexploreviz_export import metexploreviz_export

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


def extract_data_sbml(sbml_filepath):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_filepath)
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

    return genes, id_compounds, id_reactions


def test_pgdb_to_padmet():

    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    assert test_padmetSpec.dicOfNode['FAO-PWY'].misc['COMMON-NAME'][0] == 'fatty acid &beta;-oxidation I - (generic)'

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))


def test_pgdb_to_padmet_cli():

    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    test_padmetSpec = PadmetSpec('test.padmet')

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    assert test_padmetSpec.dicOfNode['FAO-PWY'].misc['COMMON-NAME'][0] == 'fatty acid &beta;-oxidation I - (generic)'

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))

    os.remove('test.padmet')


def test_pgdb_to_padmet_without_genes():
    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb')
    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)
    assert all_pwys == FABO_PWYS

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert all_genes == []


def test_pgdb_to_padmet_without_genes_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet'])

    test_padmetSpec = PadmetSpec('test.padmet')

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)
    assert all_pwys == FABO_PWYS

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert all_genes == []

    os.remove('test.padmet')


def test_pgdb_to_padmet_no_orphan_with_genes():
    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', no_orphan=True, extract_gene=True)

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    # 2.3.1.49-RXN is added as it has no association with genes
    assert '2.3.1.49-RXN' not in all_rxns

    assert set(FABO_CPDS).issubset(set(all_cpds))

    # b24589 is a manually added genes without association to a reaction
    assert 'b24589' not in all_genes


def test_pgdb_to_padmet_no_orphan_with_genes_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene', '--no-orphan'])

    test_padmetSpec = PadmetSpec('test.padmet')

    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(test_padmetSpec)

    assert all_pwys == FABO_PWYS

    # 2.3.1.49-RXN is added as it has no association with genes
    assert '2.3.1.49-RXN' not in all_rxns

    assert set(FABO_CPDS).issubset(set(all_cpds))

    # b24589 is a manually added genes without association to a reaction
    assert 'b24589' not in all_genes

    os.remove('test.padmet')


def test_pgdb_to_padmet_prot_ids70_cli():
    # Get data
    subprocess.call(
        ['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--prot-ids70'])
    test_padmetSpec = PadmetSpec('test.padmet')
    all_xref = dict()
    for node in test_padmetSpec.dicOfNode.values():
        if node.type == "xref":
            all_xref[node.id] = node.misc

    # Test test.padmet content
    assert set(all_xref['ACYLCOADEHYDROG-RXN_xrefs']['UNIPROT_70']) == {'Q47146'}
    assert set(all_xref['ENOYL-COA-DELTA-ISOM-RXN_xrefs']['UNIPROT_70']) == {'O75521', 'Q6NYL3', 'P42126', 'Q489W3',
                                                                             'P77399', 'P55100', 'O23299', 'B5XYH0',
                                                                             'O04469', 'Q8W1L6', 'Q9ZCZ1'}
    assert set(all_xref['ENOYL-COA-DELTA-ISOM-RXN_xrefs']['PID_70']) == {'AAC83700.1'}

    # Test proteins_seq_ids_reduced_70.fasta file
    files_records = list(SeqIO.parse("proteins_seq_ids_reduced_70.fasta", "fasta"))
    assert len([record.id for record in files_records]) == 24
    assert files_records[0].id == 'UNIPROT:Q58DM8'
    assert files_records[0].seq == 'MAALRALLPRVRAPLRPWLFCPVQRSFASSAAFEYIITAKKGRNSNVGLIQLNRPKALNALCNGLIVELNQALQAFEE' \
                                   'DPAVGAIVLTGGEKVFAAGADIKEMQSLTFQNCYSGGFLSHWDQLTRVKKPVIAAVNGYALGGGCELAMMCDIIYAGE' \
                                   'KAQFGQPEILIGTIPGAGGTQRLTRAVGKSLAMEMVLTGDRISAQDAKQAGLVSKIFPVETVVEEAIQCAEKIASNSK' \
                                   'IVTAMAKESVNAAFEMTLAEGVKLEKKLFYSTFATEDRKEGMAAFVEKRKANFKDQ'
    assert files_records[-1].id == 'UNIPROT:Q9P4U7'

    # Delete created files
    os.remove('test.padmet')
    os.remove('proteins_seq_ids_reduced_70.fasta')


def test_sbmlGenerator():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')

    genes, id_compounds, id_reactions = extract_data_sbml('fabo.sbml')

    assert set(FABO_RXNS).issubset(set(id_reactions))

    assert set(FABO_CPDS).issubset(set(id_compounds))

    assert set(FABO_GENES).issubset(set(genes))

    os.remove('fabo.sbml')


def test_sbmlGenerator_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])
    fabo_padmetSpec = PadmetSpec('test.padmet')

    subprocess.call(['padmet', 'sbmlGenerator', '--padmet', 'test.padmet', '--output', 'fabo.sbml', '--sbml_lvl', '3'])

    genes, id_compounds, id_reactions = extract_data_sbml('fabo.sbml')

    assert set(FABO_RXNS).issubset(set(id_reactions))

    assert set(FABO_CPDS).issubset(set(id_compounds))

    assert set(FABO_GENES).issubset(set(genes))

    os.remove('test.padmet')
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


def test_sbml_to_padmet_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'sbmlGenerator', '--padmet', 'test.padmet', '--output', 'fabo.sbml', '--sbml_lvl', '3'])

    subprocess.call(['padmet', 'sbml_to_padmet', '--sbml', 'fabo.sbml', '--padmetSpec', 'fabo.padmet'])

    fabo_padmet = PadmetSpec('fabo.padmet')
    all_pwys, all_cpds, all_rxns, all_genes = extract_data_padmet(fabo_padmet)

    assert all_pwys == []

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))

    os.remove('fabo.sbml')
    os.remove('test.padmet')
    os.remove('fabo.padmet')


def test_sbml_to_sbml():
    fabo_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(fabo_padmetSpec, 'fabo.sbml')
    from_sbml_to_sbml('fabo.sbml', 'fabo_2.sbml', 2, cpu=1)

    sbml_3_genes, sbml_3_compounds, sbml_3_reactions = extract_data_sbml('fabo.sbml')

    sbml_2_genes, sbml_2_compounds, sbml_2_reactions = extract_data_sbml('fabo_2.sbml')

    assert set(FABO_RXNS).issubset(set(sbml_2_reactions))

    assert set(FABO_CPDS).issubset(set(sbml_2_compounds))

    assert set(FABO_GENES).issubset(set(sbml_2_genes))

    assert set(sbml_3_reactions).issubset(set(sbml_2_reactions))

    assert set(sbml_3_compounds).issubset(set(sbml_2_compounds))

    assert set(sbml_3_genes).issubset(set(sbml_2_genes))

    os.remove('fabo.sbml')
    os.remove('fabo_2.sbml')


def test_sbml_to_sbml_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'sbmlGenerator', '--padmet', 'test.padmet', '--output', 'fabo.sbml', '--sbml_lvl', '3'])

    subprocess.call(['padmet', 'sbml_to_sbml', '--input', 'fabo.sbml', '--output', 'fabo_2.sbml', '--new_sbml_lvl', '2', '--cpu', '1'])

    sbml_3_genes, sbml_3_compounds, sbml_3_reactions = extract_data_sbml('fabo.sbml')

    sbml_2_genes, sbml_2_compounds, sbml_2_reactions = extract_data_sbml('fabo_2.sbml')

    assert set(FABO_RXNS).issubset(set(sbml_2_reactions))

    assert set(FABO_CPDS).issubset(set(sbml_2_compounds))

    assert set(FABO_GENES).issubset(set(sbml_2_genes))

    assert set(sbml_3_reactions).issubset(set(sbml_2_reactions))

    assert set(sbml_3_compounds).issubset(set(sbml_2_compounds))

    assert set(sbml_3_genes).issubset(set(sbml_2_genes))

    os.remove('fabo.sbml')
    os.remove('fabo_2.sbml')
    os.remove('test.padmet')


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


def test_wikiGenerator_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'fabo.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'wikiGenerator', '--padmet', 'fabo.padmet', '--output', 'output', '--wiki_id', 'TEST'])

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


def test_sbml_to_curation_form_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'sbmlGenerator', '--padmet', 'test.padmet', '--output', 'fabo.sbml', '--sbml_lvl', '3'])

    rxns = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN']
    id_reactions = ['R_'+sbmlPlugin.convert_to_coded_id(reaction) for reaction in rxns]

    with open('reactions.txt', 'w') as tmp_file:
        for id_reaction in id_reactions:
            tmp_file.write(id_reaction+'\n')

    subprocess.call(['padmet', 'sbml_to_curation_form', '--sbml', 'fabo.sbml', '--output', 'form.txt', '--rxn_file', 'reactions.txt'])

    os.remove('test.padmet')
    os.remove('fabo.sbml')
    os.remove('reactions.txt')

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


def test_gbk_to_fasta_cli():
    subprocess.call(['padmet', 'gbk_to_faa', '--gbk', 'test_data/gbk/fatty_acid_beta_oxydation_I_1.gbk',
                        '--output', 'fatty_acid_beta_oxydation_I_1.faa', '--qualifier', 'locus_tag'])

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


def test_gene_to_targets_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'gene_to_targets', '--padmetSpec', 'test.padmet', '--genes', 'test_data/genes_to_targets.txt', '--output', 'targets.txt'])

    os.remove('test.padmet')

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


def test_padmet_to_padmet_cli():
    # Using inpu data, create 2 padmets and delete one reaction in each.
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])
    fabo_1_padmetSpec = PadmetSpec('test.padmet')
    fabo_1_padmetSpec.delNode('ACYLCOASYN-RXN')
    fabo_1_padmetSpec.generateFile('fabo_1.padmet')

    _, _, all_rxns, _ = extract_data_padmet(fabo_1_padmetSpec)

    os.remove('test.padmet')

    assert not set(FABO_RXNS).issubset(set(all_rxns))

    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])
    fabo_2_padmetSpec = PadmetSpec('test.padmet')
    fabo_2_padmetSpec.delNode('ACYLCOADEHYDROG-RXN')
    fabo_2_padmetSpec.generateFile('fabo_2.padmet')

    _, _, all_rxns, _ = extract_data_padmet(fabo_2_padmetSpec)

    os.remove('test.padmet')

    assert not set(FABO_RXNS).issubset(set(all_rxns))

    # By merging them we should retrieve the two deleted reactions
    padmet_to_padmet('fabo_1.padmet,fabo_2.padmet', 'fabo.padmet')
    subprocess.call(['padmet', 'padmet_to_padmet', '--to_add', 'fabo_1.padmet,fabo_2.padmet', '--output', 'fabo.padmet'])

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


def test_padmet_to_padmet_reversibility_cli():
    """
    Test an issue encountered when a reaction is defined in a direction in one padmet.
    And in another the same reaction is defined as reversible.
    In old padmet this leads to the reaction having all of its reactants/products as reactants and also as products.
    """
    # Read padmet file
    subprocess.call(['padmet', 'padmet_to_padmet', '--to_add', 'test_data/padmet', '--output', 'fabo.padmet'])

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


def test_padmet_to_matrix_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet'])

    subprocess.call(['padmet', 'padmet_to_matrix', '--padmet', 'test.padmet', '--output', 'matrix.tsv'])

    os.remove('test.padmet')

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


def test_extract_rxn_with_gene_assoc_cli():
    subprocess.call(['padmet', 'pgdb_to_padmet', '--pgdb', 'test_data/pgdb', '--output', 'test.padmet', '--extract-gene'])

    subprocess.call(['padmet', 'sbmlGenerator', '--padmet', 'test.padmet', '--output', 'fabo.sbml'])

    # Extract reactions with only genes association so 2.3.1.49-RXN should not be here.
    subprocess.call(['padmet', 'extract_rxn_with_gene_assoc', '--sbml', 'fabo.sbml', '--output', 'fabo_rxn_with_genes.sbml'])


    reader = libsbml.SBMLReader()
    document = reader.readSBML('fabo_rxn_with_genes.sbml')
    model = document.getModel()
    reactions = model.getListOfReactions()
    id_reactions = [sbmlPlugin.convert_from_coded_id(reaction.id)[0] for reaction in reactions]

    os.remove('test.padmet')
    os.remove('fabo.sbml')
    os.remove('fabo_rxn_with_genes.sbml')

    assert '2.3.1.49-RXN' not in id_reactions

