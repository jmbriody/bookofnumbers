import pytest
from collections import namedtuple
from bookofnumbers import *


def test_quin():
    # a = qmc(2078)
    Term = namedtuple('Term', 'termset used ones source generation final')
    assert quinemc(2078) == "A'BC'D' + A'B'C + A'B'D + B'CD"
    assert quinemc(2077) == "A'C'D' + A'B'D' + B'CD"
    assert quinemc(12309, True) == "A'C'D' + A'B'D' + ABC'"
    assert quinemc(2003) == "A'C'D' + AB'D' + A'BC + B'C'"
    assert quinemc(255) == "1"
    assert quinemc(0) == "0"

    canon_list = ["ABC'D", "A'B'CD'", "ABC'D'", "A'BC'D'", "A'B'C'D'"]
    canon_string = "ABC'D + A'B'CD' + ABC'D' + A'BC'D' + A'B'C'D'"

    canon_string_error = "ABCD + A'B'D' + ABC'D' + A'BC'D' + A'B'C'D'"
    assert quinemc(canon_string) == "A'B'D' + A'C'D' + ABC'"
    assert quinemc(canon_list) == "A'B'D' + A'C'D' + ABC'"

    assert isinstance(quinemc(canon_string_error), ValueError)
    
    second_result = ["AD'", "B'C'D", "A'BC'"]
    a, b, c = quinemc(2046, False, True)
    c1 = c[1]
    r = ["".join(sorted(ti.termset)) for ti in c1]
    assert set(r) == set(second_result)
    assert len(b) == 26

