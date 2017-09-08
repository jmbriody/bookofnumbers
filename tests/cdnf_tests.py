import pytest
from bookofnumbers import *

def test_quin():
    # a = qmc(2078)
    assert quinmc(2078) == "A'BC'D' + A'B'C + A'B'D + B'CD"
    assert quinmc(2077) == "A'C'D' + B'CD + A'B'D'"
    assert quinmc(12309, True) == "A'B'D' + ABC' + A'C'D'"
    assert quinmc(2003) == "AB'D' + A'BC + A'C'D' + B'C'"
    assert quinmc(255) == "1"
    assert quinmc(0) == "0"

    canon_list = ["ABC'D", "A'B'CD'", "ABC'D'", "A'BC'D'", "A'B'C'D'"]
    canon_string = "ABC'D + A'B'CD' + ABC'D' + A'BC'D' + A'B'C'D'"

    canon_string_error = "ABCD + A'B'D' + ABC'D' + A'BC'D' + A'B'C'D'"
    assert quinmc(canon_string) == "A'B'D' + ABC' + A'C'D'"
    assert quinmc(canon_list) == "A'B'D' + ABC' + A'C'D'"

    assert isinstance(quinmc(canon_string_error), ValueError)
    



