from padmet.utils import gbr

def test_gbr():
    expected = [('a', 'b', 'd'), ('a', 'b', 'e'), ('a', 'c', 'd'), ('a', 'c', 'e')]
    gbr_results = [elements for elements in gbr.compile_input('a&(b|c)&(d|e)')]

    assert gbr_results == expected
