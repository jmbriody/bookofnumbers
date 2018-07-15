# CDNF

[![Build Status](https://travis-ci.org/jmbriody/cdnf.svg?branch=master)](https://travis-ci.org/jmbriody/cdnf)
[![Coverage Status](https://coveralls.io/repos/github/jmbriody/cdnf/badge.svg?branch=master)](https://coveralls.io/github/jmbriody/cdnf?branch=master)

# CDNF and Quine McCluskey in Python
So, this started as a simple app to take an int and create its relevant CDNF (Canonical Disjunctive Noraml Form) based on Section 3.4 of "Digital Design" from Randall Hyde's "The Art of Assembly Language Programming". When that proved to be too easy I added a quine mccluskey implementation. Unlike other implementations I've seen results are returned using Alpha characters. E.g. "ABD' + C" rather than "11-0 + --1-". Sorry.

## Installation
There are no requirements and no extra stuff needed--it is pure python. Should work with 2.7, 3.4, 3.5, and 3.6. To install you can . . .

git clone git@github.com:jmbriody/cdnf.git

cd cdnf

python setup.py

... OR ...

pip install git+git@github.com:jmbriody/cdnf.git

. . . As always--do this crud in a virtual environment. 

To use just do `from cdnf import *`

## The things it does
canonical(item, highorder_a, includef):

    Where it all began. Takes and int (item) and returns the Canonical Disjunctive Normal Form for that int. So using 248 as an example (binary version is 1111000) the CDNF is "ABC' + ABC + AB'C' + AB'C + A'BC". Getting that is `canonical(248)`.

    highorder_a -- determines whether A is the high order bit (default) or low order bit.

    includef -- determines whether the input is returned as part of the result. 

to_cdnf(item, ranged=False):

    This is where we get funky. This will take a minimized function and return the CDNF version of it. Using "ranged" will add missing characters. So . . .
    to_cdnf("AB + D") will return "ABD' + ABD + AB'D + A'BD + A'B'D" (note no "C").
    to_cdnf("AB + D", ranged=True) will return  "ABCD' + ABCD + ABC'D' + ABC'D + AB'CD + AB'C'D + A'BCD + A'BC'D + A'B'CD + A'B'C'D"

    Using ranged can be dangerous `to_cdnf("A + I", ranged=True)` returns 384 terms; `to_cdnf("A + Z", ranged=True)` returns 50,000,000 terms. 

quinemc(myitem, highorder_a=True, full_results=False):

    Short for Quinne-McCluskey algorithm:
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms. Terms can be
    composed of any alpha characters (plus (')).

    highorder_a: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True three items are returned:
    --result: Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible results. A sparse few set of equations can be reduced to
    other equations of the same length. This is a dictionary of those possible alternatives for the
    occasional items where this occurs. For example with quinemc(743) all of ...
    "B'C'D + A'B'D' + A'CD' + A'BD"
    "B'C'D + A'B'D' + A'C'D + A'BC",
    "B'C'D + A'B'D' + A'BC + A'BD",
    "B'C'D + A'B'C' + A'CD' + A'BD"
    are equivalent reductions.

## Helper functions
result_to_int(t):

    If you run `r, t, p = quinemc(to_cdnf("B'C'D + A'CD' + A'BD + A'B'D'"), full_results=True)`
    result_to_int(t) will equal 743

alternatives(t, p):

    Similarly, running `r, t, p = quinemc(to_cdnf("B'C'D + A'CD' + A'BD + A'B'D'"),
    full_results=True)`, then alternatives(t, p) will provide a list of equivalent reduced
    forms. (Majority of items do **not** have alternatives 743 just happens to have some).
    In this case (for 743) the main result (r) is "B'C'D + A'CD' + A'BD + A'B'D'"
    alternatives(t,p) is . . .
    ["B'C'D + A'B'D' + A'CD' + A'BD",
     "B'C'D + A'B'D' + A'C'D + A'BC",
     "B'C'D + A'B'D' + A'BC + A'BD",
     "B'C'D + A'B'C' + A'CD' + A'BD"]


## An interesting example
r, s, t = quinemc(743, full_results=True)

r will be: "B'C'D + A'CD' + A'BD + A'B'D'"

s . . .

[Term(termset={"B'", "A'", "D'", "C'"}, used=True, ones=0, source=[0], generation=1, final=None, binary='0000', row=0, dontcare=None),
 Term(termset={"B'", "A'", "D'", 'C'}, used=True, ones=1, source=[1], generation=1, final=None, binary='0010', row=1, dontcare=None),
 Term(termset={"B'", "A'", 'D', "C'"}, used=True, ones=1, source=[2], generation=1, final=None, binary='0001', row=2, dontcare=None),
 Term(termset={"B'", 'D', "C'", 'A'}, used=True, ones=2, source=[3], generation=1, final=None, binary='1001', row=3, dontcare=None),
 Term(termset={'C', "A'", "D'", 'B'}, used=True, ones=2, source=[4], generation=1, final=None, binary='0110', row=4, dontcare=None),
 Term(termset={"A'", 'D', "C'", 'B'}, used=True, ones=2, source=[5], generation=1, final=None, binary='0101', row=5, dontcare=None),
 Term(termset={'C', "A'", 'D', 'B'}, used=True, ones=3, source=[6], generation=1, final=None, binary='0111', row=6, dontcare=None),
 Term(termset={"B'", "A'", "D'"}, used=False, ones=0, source=[0, 1], generation=2, final='Added', binary=None, row=7, dontcare=None),
 Term(termset={"B'", "A'", "C'"}, used=False, ones=0, source=[0, 2], generation=2, final=None, binary=None, row=8, dontcare=None),
 Term(termset={'C', "A'", "D'"}, used=False, ones=1, source=[1, 4], generation=2, final='Added', binary=None, row=9, dontcare=None),
 Term(termset={"B'", 'D', "C'"}, used=False, ones=1, source=[2, 3], generation=2, final='Required', binary=None, row=10, dontcare=None),
 Term(termset={"A'", 'D', "C'"}, used=False, ones=1, source=[2, 5], generation=2, final=None, binary=None, row=11, dontcare=None),
 Term(termset={'C', "A'", 'B'}, used=False, ones=2, source=[4, 6], generation=2, final=None, binary=None, row=12, dontcare=None),
 Term(termset={"A'", 'D', 'B'}, used=False, ones=2, source=[5, 6], generation=2, final='Added', binary=None, row=13, dontcare=None)]

t . . .

defaultdict(list,
            {0: [Term(termset={"B'", "A'", "D'"}, used=False, ones=0, source=[0, 1], generation=2, final='Added', binary=None, row=7, dontcare=None),
              Term(termset={'C', "A'", "D'"}, used=False, ones=1, source=[1, 4], generation=2, final='Added', binary=None, row=9, dontcare=None),
              Term(termset={"A'", 'D', 'B'}, used=False, ones=2, source=[5, 6], generation=2, final='Added', binary=None, row=13, dontcare=None)],
             1: [Term(termset={"B'", "A'", "D'"}, used=False, ones=0, source=[0, 1], generation=2, final='Added', binary=None, row=7, dontcare=None),
              Term(termset={"A'", 'D', "C'"}, used=False, ones=1, source=[2, 5], generation=2, final=None, binary=None, row=11, dontcare=None),
              Term(termset={'C', "A'", 'B'}, used=False, ones=2, source=[4, 6], generation=2, final=None, binary=None, row=12, dontcare=None)],
             2: [Term(termset={"B'", "A'", "D'"}, used=False, ones=0, source=[0, 1], generation=2, final='Added', binary=None, row=7, dontcare=None),
              Term(termset={'C', "A'", 'B'}, used=False, ones=2, source=[4, 6], generation=2, final=None, binary=None, row=12, dontcare=None),
              Term(termset={"A'", 'D', 'B'}, used=False, ones=2, source=[5, 6], generation=2, final='Added', binary=None, row=13, dontcare=None)],
             3: [Term(termset={"B'", "A'", "C'"}, used=False, ones=0, source=[0, 2], generation=2, final=None, binary=None, row=8, dontcare=None),
              Term(termset={'C', "A'", "D'"}, used=False, ones=1, source=[1, 4], generation=2, final='Added', binary=None, row=9, dontcare=None),
              Term(termset={"A'", 'D', 'B'}, used=False, ones=2, source=[5, 6], generation=2, final='Added', binary=None, row=13, dontcare=None)]})



