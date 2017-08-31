import pytest
from bookofnumbers.cdnf import quinmc

def test_quin():
    # a = qmc(2078)
    assert quinmc(2078) == "A'B'CD' + BC'D' + AC'D' + ABC'"
    assert quinmc(2077) == "A'B'D' + ABC' + A'C'D'"
    assert quinmc(2003) == "A'C'D + BCD' + A'B'D' + B'C'"
    assert quinmc(255) == "1"
    assert quinmc(0) == "0"

    canon_list = ["ABC'D", "A'B'CD'", "ABC'D'", "A'BC'D'", "A'B'C'D'"]
    canon_string = "ABC'D + A'B'CD' + ABC'D' + A'BC'D' + A'B'C'D'"

    canon_string_error = "ABCD + A'B'D' + ABC'D' + A'BC'D' + A'B'C'D'"
    assert quinmc(canon_string) == "A'B'D' + ABC' + A'C'D'"
    assert quinmc(canon_list) == "A'B'D' + ABC' + A'C'D'"

    assert isinstance(quinmc(canon_string_error), ValueError)
    



