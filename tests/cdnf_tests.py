import pytest
from collections import namedtuple
from cdnf import *


def test_quin():
    # a = qmc(2078)
    assert isinstance(canonical("ABC"), ValueError)
    assert canonical(2077, True, True) == "f(2077) = AB'CD + A'BC'D' + A'B'CD' + A'B'CD + A'B'C'D'"
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
    assert isinstance(quinemc({"ABC", "A'C"}), ValueError)

    assert to_cdnf("B'CD + A'C'D' + A'B'D'") == "AB'CD + A'BC'D' + A'B'CD' + A'B'CD + A'B'C'D'"
    assert to_cdnf(["B'CD", "A'C'D'", "A'B'D'"]) == "AB'CD + A'BC'D' + A'B'CD' + A'B'CD + A'B'C'D'"
    assert to_cdnf("C + A") == "AC' + AC + A'C"
    assert to_cdnf("C + A", 1) == "ABC' + ABC + AB'C' + AB'C + A'BC + A'B'C"
    assert to_cdnf("ry + t") == "rty' + rty + rt'y + r'ty' + r'ty"

    assert isinstance(to_cdnf(2077), ValueError)

    assert isinstance(quinemc(canon_string_error), ValueError)
    
    second_result = ["A'D", "AB'C'", "B'CD'"]
    a, b, c = quinemc(2046, True, True)
    c1 = c[1]
    r = ["".join(sorted(ti.termset)) for ti in c1]
    assert set(r) == set(second_result)
    assert len(b) == 26

    r, s, t = quinemc(743, 1, 1)
    assert result_to_int(s) == 743
    assert alternatives(s, t)[2] == "B'C'D + A'B'D' + A'BC + A'BD"

    assert quinemc([743, [0, 1]]) == "B'C'D + A'CD' + A'BD"




