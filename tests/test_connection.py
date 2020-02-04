import libsbml
import os

from padmet.utils import sbmlPlugin
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet.utils.connection.sbmlGenerator import padmet_to_sbml

FABO_RXNS = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN',
            'ENOYL-COA-DELTA-ISOM-RXN', 'OHBUTYRYL-COA-EPIM-RXN', 'KETOACYLCOATHIOL-RXN',
            'OHACYL-COA-DEHYDROG-RXN']

FABO_PWYS = ['FAO-PWY']

FABO_CPDS = ['D-3-HYDROXYACYL-COA', 'L-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA',
            'CIS-DELTA3-ENOYL-COA', 'Saturated-Fatty-Acyl-CoA', 'ETF-Oxidized',
            'ETF-Reduced', 'WATER', 'CPD0-1163', 'CPD0-1162', 'CPD0-1158',
            'PROTON', 'CPD66-39', 'CO-A', 'ATP', 'PPI', 'AMP', 'ACETYL-COA',
            '3-KETOACYL-COA', 'NAD', 'NADH']

FABO_GENES = ['b2341', 'b0221', 'b1805', 'b1701', 'b2342', 'b3845', 'b3846']


def test_pgdb_to_padmet():

    padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    total_pwy_id = set()
    total_cpd_id = set()

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node for node in padmetSpec.dicOfNode.values() if node.type == "gene"]
    nb_rxn_with_ga = 0
    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        if any([rlt for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_linked_to"]):
            nb_rxn_with_ga += 1
        total_pwy_id.update(pathways_ids)

    all_pwys = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node_id for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]
    all_rxns = [node.id for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node.id for node in padmetSpec.dicOfNode.values() if node.type == "gene"]


    assert all_pwys == FABO_PWYS

    assert set(FABO_RXNS).issubset(set(all_rxns))

    assert set(FABO_CPDS).issubset(set(all_cpds))

    assert set(FABO_GENES).issubset(set(all_genes))


def test_sbmlGenerator():
    padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)
    padmet_to_sbml(padmetSpec, 'fabo.sbml')
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