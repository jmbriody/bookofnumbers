"""Functions for generating boolean algebra normal forms and minimization

Provides two features. 1) Takes an integer and returns the Canonical Normal Form associated
with that integer (based on Section 3.4 of "Digital Design" from Randall Hyde's "The Art of
Assembly Language Programming"). And 2) Will take a canonical logic statement and return a
minimized version.

General format of results are:
    "AB" means A AND B
    "C'" is Not C
    "+" is Or
So "AB + C'" is (A AND B) or (Not C)

In the standard mode (highorder_a = True) A is the high order bit:

    BP  |  A  |  B  |  C     248 (11111000)
    ---------------------    ---
    1   |  0  |  0  |  0      0
    2   |  0  |  0  |  1      0
    3   |  0  |  1  |  0      0
    4   |  0  |  1  |  1      1     A'BC
    5   |  1  |  0  |  0      1     AB'C'
    6   |  1  |  0  |  1      1     AB'C
    7   |  1  |  1  |  0      1     ABC'
    8   |  1  |  1  |  1      1     ABC

When highorder_a is False A would be the low order bit for 248.

Two primary user functions are:
caononical(x, highorder_a=False, includef=False):
    Simplest call is canonical(248) which will return the canonical boolean algebraic expression
    for "248". In the default setup A is the Low Order bit (highorder_a=False)--i.e. 0101.
    When highorder_a is True A will be the High Order bit--i.e. 0000111100001111. includef
    determines whether "F(248) = " should be prepended to the result.



>>> canonical(248)
"ABC + ABC' + AB'C + AB'C' + A'BC"
>>> canonical(248, False)
"ABC + A'BC + AB'C + A'B'C + ABC'"
>>> canonical(248, True, True)
"f(248) = ABC + ABC' + AB'C + AB'C' + A'BC"

quinemc(myitem, highorder_a=True, full_results=False):
    Short for Quinn-McCluskey algorithm:
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms

    highorder_a: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --final_result: List of Named Tuples "termset" containing terms of the Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method. May be empty.

>>> quinemc(248)
'BC + A'
>>> quinemc(248, False)
'AB + C'
>>> quinemc("ABC + A'BC + AB'C + A'B'C + ABC'")
'AB + C'
>>> quinemc(["ABC", "A'BC", "AB'C", "A'B'C", "ABC'"])
'AB + C'
>>> result, final_result, term_list, possibles = quinemc(248, False, True)

"""
from __future__ import division  # needed for python2
import re
import math
import itertools
import string
from collections import namedtuple, defaultdict
from operator import attrgetter
from profilehooks import coverage, profile

# Term is a namedtuple used by the Quin-McCluskey reduction portion of the code.
# A list of Terms is used for the minimize process and another list is used for
# the final list of minimized terms
#
# --termset: A set of charcaters (e.g. ("A'", 'B', 'D')) representing a minterm
# --used: True/False depending on whether this term has been used in a later generation
# --ones: Number of un-primed characters. Used for the reduction process.
# --source: list of indexes of terms used to derive the current term
# --generation: which generation the term belongs to
#
# For minimize/reduction on items in the latest generation are compared. With in a
# generation only terms where "ones" differs by 1 can be merged/minimized.
#
# Code relies on doing set comparisons of characters as opposed to using 1's and 0's
# and position holders (e.g. AB'D vs. 10-1). For very large reductions (e.g.
# quinemc(4222345678921334)) the number of permutations and comparisons grows pretty
# large and can be slow.
Term = namedtuple('Term', 'termset used ones source generation final')

@coverage
def canonical(x, highorder_a=True, includef=False):
    """
    Takes an integer and returns the boolean algebra canonical disjunctive normal form.

    Arguments:
    x: integer to convert to a canonical boolean expression. based on Section 3.4 of
    "Digital Design" from Randall Hyde's "The Art of Assembly Language Programming"
http://www.plantation-productions.com/Webster/www.artofasm.com/Linux/PDFs/DigitalDesign.pdf

    highorder_a: When "False" (default) A is the low order bit. When True A is the high order
    bit.

    includef: When True will append "F(x) = " before the result string

    """
    binary = format(x, 'b')
    binary = binary[::-1]
    letters = len(format(len(binary) - 1, 'b'))
    if letters == 1:
        letters = 2
    # highorder_a the binary and create a list of all the positions with a "1".

    indexes = [(format(i.start(), '0' + str(letters) + 'b')) for i in re.finditer('1', binary)]

    miniterms = [_minterms_(m, highorder_a) for m in indexes[::-1]]
    result = ' + '.join(miniterms)

    if includef is True:
        result = "f(" + str(x) + ") = " + result
    return result

@coverage
def _minterms_(m, highorder_a):
    alpha = sorted(string.ascii_letters)
    result = ''
    if highorder_a is False:
        m = m[::-1]

    for i, x in enumerate(m):
        result += alpha[i]
        if x == '0':
            result += "'"
    return result


# --- END OF Canonical functions ---

@coverage
def quinemc(myitem, highorder_a=True, full_results=False):
    '''Short for Quinn-McCluskey algorithm
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    Arguments:
    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms

    highorder_a: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --final_result: List of Named Tuples "termset" containing terms of the Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method. May be empty.

    >>> quinemc(248)
    'AB + C'
    >>> quinemc(248, True)
    'BC + A'
    >>> quinemc("ABC + A'BC + AB'C + A'B'C + ABC'")
    'AB + C'
    >>> quinemc(["ABC", "A'BC", "AB'C", "A'B'C", "ABC'"])
    'AB + C'
    >>> result, final_result, term_list, possibles = quinemc(248, False, True)

    Invalid input will return a ValueError (e.g. quinemc(["A'BC", "AB"]) is invalid
    because second term must contain a "C".
    '''
    if isinstance(myitem, int):
        cdnf = canonical(myitem, highorder_a).split(' + ')
    elif isinstance(myitem, str):
        cdnf = re.split(r"[^a-zA-Z']+", myitem)
    elif isinstance(myitem, list):
        cdnf = myitem
    else:
        return ValueError(myitem, "Not valid input")

    test_string = "".join(sorted(re.sub("'", '', cdnf[0])))
    for item in cdnf:
        if test_string != "".join(sorted(re.sub("'", '', item))):
            return ValueError("Term: ", item, " doesn't match valid test ", test_string)

    if full_results:
        return _minimize_(cdnf)
    else:
        return _minimize_(cdnf)[0]


@coverage
def _minimize_(cdnf):
    """
    Basic driver for performing the over all minimization. Most of the "heavy" work
    is done in _implicants_()

    """
    done = False
    current_generation = 1
    # Step 1: take terms from the caononical form and putting them into generation one
    # of our term_list
    term_list = _create_first_generation_(cdnf)

    # Step 2: merge terms of each generation to create next generation until no more merges
    # are possible (_merge_terms_ and _create_new_tuples_)
    while not done:
        done = _merge_terms_(term_list, current_generation)
        current_generation += 1

    # Step 3: Generate our final result from all terms in term_list that have not been used in
    # a merge
    possibles = _implicants_(term_list)

    result = ["".join(sorted(tempItem.termset)) for tempItem
              in term_list
              if tempItem.final is not None]
    result = " + ".join(result)

    # Handle special cases-- quinemc(0), quinemc(15), quinemc(255),
    # quinemc("AB + A'B + AB' + A'B'"), etc.
    if len(term_list) == 1 and len(term_list[0].termset) == 0:
        result = "0"
    if result == "":
        result = "1"

    return result, term_list, possibles


@coverage
def _create_first_generation_(terms):
    my_letters = set(sorted(string.ascii_letters))
    temp_terms = [set(re.findall("([A-Za-z]'*)", x)) for x in terms]  # Convert to list of sets

    # Remove duplicate terms if called with something like quinemc("ABCD + CDBA + ABC'D + DC'AB")
    temp_terms = list(temp_terms for temp_terms, _ in itertools.groupby(temp_terms))

    temp_list = [Term(x, False, len(x.intersection(my_letters)), None, 1, None)
                 for x in temp_terms]
    temp_list = sorted(temp_list, key=attrgetter('ones'))
    for idx, item in enumerate(temp_list):
        temp_list[idx] = item._replace(source=[idx])
    return temp_list


@coverage
def _merge_terms_(term_list, gen):
    done = False
    new_terms = _create_new_terms_(term_list, gen)

    if len(new_terms) > 0:
        term_list.extend(new_terms)
    else:
        done = True

    return done


@coverage
def _create_new_terms_(orig_term_list, gen):
    used_dict = {}  # a dictionary for used items
    sources = []  # avoid duplicate merges
    result = []
    # my_letters = set(string.ascii_letters)
    working_list = [x for x in orig_term_list if x.generation == gen]

    for x in working_list:
        for y in working_list:
            if y.ones == (x.ones + 1):
                sym_set = y.termset.symmetric_difference(x.termset)

                if len(sym_set) == 2 and list(sym_set)[0][0] == list(sym_set)[1][0]:
                    x_index, y_index = orig_term_list.index(x), orig_term_list.index(y)
                    used_dict[x_index] = True
                    used_dict[y_index] = True
                    new_term = y.termset.intersection(x.termset)
                    source = sorted(y.source + x.source)
                    if source not in sources:
                        sources.append(source)
                        result.append(Term(new_term, False,
                                           len(new_term.intersection(set(string.ascii_letters))),
                                           source, (gen + 1), None))
    result = sorted(result, key=attrgetter('ones'))

    for key in used_dict.keys():
        current = orig_term_list[key]
        orig_term_list[key] = Term(current.termset, True, current.ones,
                                   current.source, current.generation, None)

    return result


# ---- Finding Implicants / Final Result ---- #
# _implicants_, _required_sources_, _make_find_dict_, and _check_combinations_ do the heavy
# lifting to find the final reduced form
@coverage
def _implicants_(term_list):
    '''
    Finds terms that will cover unused cases if the essential prime implicants are not
    enough.
    '''
    possible_terms = defaultdict(list)
    list_of_sources = []

    for x in [zed for zed in term_list if zed.used is False]:
        list_of_sources += list(itertools.chain(x.source))

    required = [x for x in list_of_sources if list_of_sources.count(x) == 1]
    keep_columns = _get_columns_(term_list, required)

    # if _get_columns_ ends with nothing in keep_columns it means essential prime implicants
    # are all that is needed so we are done
    finished = bool(len(keep_columns) == 0)

    # check if single term will "cover" remaining items e.g. qmc(2077)
    if not finished:
        find_dict = _make_find_dict_(term_list, keep_columns)
        for idx, val in find_dict.items():
            if set(keep_columns) == val.sourceSet:
                term_list[idx] = term_list[idx]._replace(final="Added")
                finished = True
                break

    # if a single term doesn't cover the remaining 1st gen items start looking
    # for combinations of Terms that will fit the bill.
    if not finished:
        possible_terms = _check_combinations_(find_dict, term_list, keep_columns)

    return possible_terms


@coverage
def _get_columns_(term_list, required):
    """
    needed -- list of "Term" tuples that have used==False
    required -- integers for terms that are essential prime implicants . . .
        each required int will appear in the source list for only 1 item in needed
    """
    ignore = []
    keep = []

    for j, x in [(i, k) for i, k in enumerate(term_list) if k.used is False]:
        # Find Terms in "needed" that exist in required, add them to the final result,
        # and add that Term's sources to the "columns" we can now ignore (already covered
        # terms)
        if len((set(required) & set(x.source))) >= 1:
            term_list[j] = x._replace(final="Required")
            ignore += itertools.chain(x.source)
        # Otherwise add the sources to our list of "columns" we need to keep
        else:
            keep += itertools.chain(x.source)

    # create a list of the remaining 1st gen terms that we still need to find minterms for
    keep = list(set(keep) - set(ignore))

    return keep


@coverage
def _make_find_dict_(term_list, keep_columns):
    # Creates a dictionary referencing the remaining tuples that can potentially complete
    # the minimized form.
    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = {}
    for idx, item in [(i, k) for i, k in enumerate(term_list)
                      if k.used is False and k.final is None]:
        temp_source = set(keep_columns) & set(item.source)

        if len(temp_source) > 0:
            temp_tuple = search_tuple(temp_source, len(item.termset))
            find_dict[idx] = temp_tuple

    return find_dict


@coverage
def _check_combinations_(find_dict, term_list, keep_columns):
    match_idx = []
    matches = []
    possible_terms = defaultdict(list)
    indexes = []

    for x in range(2, (len(find_dict) + 1)):
        for items in itertools.combinations(find_dict.keys(), x):
            combined_sources = set()
            temp_count = 0
            for idx in items:
                combined_sources.update(find_dict[idx].sourceSet)
                temp_count += find_dict[idx].length
            if set(keep_columns) == combined_sources:
                matches.append(items)
                match_idx.append(temp_count)

    if len(matches) > 0:
        indexes = [matches[i] for i, v in enumerate(match_idx) if v == min(match_idx)]
        for idx, value in enumerate(indexes):
            for i in value:
                if idx == 0:
                    term_list[i] = term_list[i]._replace(final="Added")
                possible_terms[idx].append(term_list[i])

    return possible_terms


if __name__ == "__main__":
    #    print("here")
    #    a, b, c, d = quinemc(42589768824798729982179, False, True)
    #    print(a)
    #    a = quinemc(2077)
    print("2003", quinemc(2003))
    print("2077", quinemc(2077))
    print("2078", quinemc(2078))
    print("2046", quinemc(2046))
    canonical(2046)

# 65024
# 15872
# 48640
# **Situation 4 -- 38400**


# a,  b,  c,  d = quinemc(2003, False, True)
#    for i,  v in d.items():
#        print(i, v)
#    mj = canonical(14)
# skokie--6:30 -- 8:45 concert
#    mk = canonical(248)
#    mi = canonical(22256)
#    print(mj, mk, mi)

#    a, b, c, d = quinemc(2078)
#    print(a)
#    for x in b:
#        print(x)
#    print("Termlist")
#    for x in c:
#        print(x)
#    print("Dictionary")
#    for k, v, in d.items():
#        print(k, ":\t", v)
