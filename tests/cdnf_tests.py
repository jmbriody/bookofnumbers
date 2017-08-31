import pytest
from bookofnumbers.cdnf import qmc, quinmc

def test_quin():
    a = qmc(2078)
    assert a == "A'B'CD' + BC'D' + AC'D' + ABC'"
