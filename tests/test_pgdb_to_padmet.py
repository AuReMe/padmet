
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet

def test_pgdb_top_admet():
    fabo_rnxs = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN',
                'ENOYL-COA-DELTA-ISOM-RXN', 'OHBUTYRYL-COA-EPIM-RXN', 'KETOACYLCOATHIOL-RXN',
                'OHACYL-COA-DEHYDROG-RXN']

    fabo_pwy = ['FAO-PWY']


    fabo_cpds = ['D-3-HYDROXYACYL-COA', 'L-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA',
                'CIS-DELTA3-ENOYL-COA', 'Saturated-Fatty-Acyl-CoA', 'ETF-Oxidized',
                'ETF-Reduced', 'WATER', 'CPD0-1163', 'CPD0-1162', 'CPD0-1158',
                'PROTON', 'CPD66-39', 'CO-A', 'ATP', 'PPI', 'AMP', 'ACETYL-COA',
                '3-KETOACYL-COA', 'NAD', 'NADH']

    fabo_genes = ['b2341', 'b0221', 'b1805', 'b1701', 'b2342', 'b3845', 'b3846']

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


    assert all_pwys == fabo_pwy

    assert set(fabo_rnxs).issubset(set(all_rxns))

    assert set(fabo_cpds).issubset(set(all_cpds))

    assert set(fabo_genes).issubset(set(all_genes))
