
from padmet.utils import sbmlPlugin

def test_sbmlPlugin():
    # test convert_from_coded_id
    coded = "R_RXN__45__11921"
    assert sbmlPlugin.convert_from_coded_id(coded) == ("RXN-11921", "R", None)

    coded = "R_R00332_c"
    assert sbmlPlugin.convert_from_coded_id(coded) == ("R00332", "R", "c")

    coded = "R_SUCCt2_2"
    assert sbmlPlugin.convert_from_coded_id(coded) == ("SUCCt2_2", "R", None)

    coded = "S_N6_45__40_L_45_1_44_3_45_Dicarboxypropyl_41__45_L_45_lysine_c"
    assert sbmlPlugin.convert_from_coded_id(coded, pattern="_", species_tag="S") == (
        "N6-(L-1,3-Dicarboxypropyl)-L-lysine",
        "S",
        "c",
    )

    coded = "S__40_2R_41__45_2_45_Hydroxy_45_3_45__40_phosphonooxy_41__45_propanal_c"
    assert sbmlPlugin.convert_from_coded_id(coded, pattern="_", species_tag="S") == (
        "(2R)-2-Hydroxy-3-(phosphonooxy)-propanal",
        "S",
        "c",
    )

    coded = "M_citr_L_m"
    assert sbmlPlugin.convert_from_coded_id(coded) == ("citr_L", "M", "m")
