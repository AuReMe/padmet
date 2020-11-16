from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet

FABO_RXNS = ['ACYLCOADEHYDROG-RXN', 'ACYLCOASYN-RXN', 'ENOYL-COA-HYDRAT-RXN',
            'ENOYL-COA-DELTA-ISOM-RXN', 'OHBUTYRYL-COA-EPIM-RXN', 'KETOACYLCOATHIOL-RXN',
            'OHACYL-COA-DEHYDROG-RXN']

FABO_PWYS = ['FAO-PWY']

FABO_CPDS = ['D-3-HYDROXYACYL-COA', 'L-3-HYDROXYACYL-COA', 'TRANS-D2-ENOYL-COA',
            'CIS-DELTA3-ENOYL-COA', 'Saturated-Fatty-Acyl-CoA', 'ETF-Oxidized',
            'ETF-Reduced', 'WATER', 'PROTON', 'CPD66-39', 'CO-A', 'ATP', 'PPI', 'AMP', 'ACETYL-COA',
            '3-KETOACYL-COA', 'NAD', 'NADH']

FABO_GENES = ['b2341', 'b0221', 'b1805', 'b1701', 'b2342', 'b3845', 'b3846']

def test_getter():
    test_padmetSpec = from_pgdb_to_padmet('test_data/pgdb', extract_gene=True)

    assert FABO_PWYS == test_padmetSpec.getPathways()

    assert set(FABO_RXNS).issubset(set(test_padmetSpec.getReactions()))

    assert set(FABO_CPDS).issubset(set(test_padmetSpec.getCompounds()))

    assert set(FABO_GENES).issubset(set(test_padmetSpec.getGenes()))
