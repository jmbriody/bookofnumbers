import pytest
from collections import namedtuple
from bookofnumbers import *


def test_quin():
    # a = qmc(2078)
    Term = namedtuple('Term', 'termset used ones source generation final')
    assert quinemc(2078) == "B'CD + A'BC'D' + A'B'D + A'B'C"
    assert quinemc(2077) == "B'CD + A'C'D' + A'B'D'"
    assert quinemc(12309, True) == "ABC' + A'C'D' + A'B'D'"
    assert quinemc(2003) == "B'C' + AB'D' + A'C'D' + A'BC"
    assert quinemc(255) == "1"
    assert quinemc(0) == "0"

    canon_list = ["ABC'D", "A'B'CD'", "ABC'D'", "A'BC'D'", "A'B'C'D'"]
    canon_string = "ABC'D + A'B'CD' + ABC'D' + A'BC'D' + A'B'C'D'"

    canon_string_error = "ABCD + A'B'D' + ABC'D' + A'BC'D' + A'B'C'D'"
    assert quinemc(canon_string) == "ABC' + A'C'D' + A'B'D'"
    assert quinemc(canon_list) == "ABC' + A'C'D' + A'B'D'"

    assert isinstance(quinemc(canon_string_error), ValueError)
    
    second_result = ["AD'", "B'C'D", "A'BC'"]
    a, b, c = quinemc(2046, False, True)
    c1 = c[1]
    r = ["".join(sorted(ti.termset)) for ti in c1]
    assert set(r) == set(second_result)
    assert len(b) == 26

