"""Functions for generating boolean algebra normal forms and minimization

Provides three features.
1) Takes an integer and returns the Canonical Normal Form associated with that
integer (based on Section 3.4 of "Digital Design" from Randall Hyde's "The Art of
Assembly Language Programming").

2) Will take a minimized equation and return "a" Canonical Normal Form

3) Will take a canonical logic statement and return a minimized version.


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

Three primary user functions are:

caononical(x, highorder_a=False, includef=False):
    Simplest call is canonical(248) which will return the canonical boolean algebraic expression
    for "248". In the default setup A is the High Order bit (highorder_a=True)--i.e. 00001111.
    When highorder_a is False A will be the Low Order bit--i.e. 01010101. includef
    determines whether "F(248) = " should be prepended to the result.


>>> canonical(248)
"ABC' + ABC + AB'C' + AB'C + A'BC"
>>> canonical(248, False)
"ABC' + ABC + AB'C + A'BC + A'B'C"
>>> canonical(248, True, True)
"f(248) = ABC' + ABC + AB'C' + AB'C + A'BC"

to_cdnf(x):


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
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method that would create a
    smaller form--items may be "longer" than the resulting minimized form.

>>> quinemc(248)
'BC + A'
>>> quinemc(248, False)
'C + AB'
>>> quinemc("ABC + A'BC + AB'C + A'B'C + ABC'")
'C + AB'
>>> quinemc(["ABC", "A'BC", "AB'C", "A'B'C", "ABC'"])
'C + AB'
>>> quinemc("rts + r'ts + rt's + r't's + rts'")
's + rt'
>>> result, term_list, possibles = quinemc(248, False, True)

"""
import re
import itertools
import string
# import pprint as pp
# import line_profiler
from collections import namedtuple, defaultdict
from operator import attrgetter, itemgetter
#from profilehooks import coverage, profile

# Term is a namedtuple used by the Quin-McCluskey reduction portion of the code.
# A list of Terms is used for the minimize process and another list is used for
# the final list of minimized terms
#
# --termset: A set of charcaters (e.g. ("A'", 'B', 'D')) representing a minterm
# --used: True/False depending on whether this term has been used in a later generation
# --ones: Number of un-primed characters. Used for the reduction process.
# --source: list of indexes of terms used to derive the current term
# --generation: which generation the term belongs to
# --final: (None|Required|Added). "None"--default setting; term isn't used in reduced
#       form. "Required"--terms that are required for the reduced form. "Added"--
#       terms that complete the reduced form using Petrick's method.
#
# For minimize/reduction items in the latest generation are compared. With in a
# generation only terms where "ones" differs by 1 can be merged/minimized.
#
# Code relies on doing set comparisons of characters as opposed to using 1's and 0's
# and position holders (e.g. AB'D vs. 10-1). For very large reductions (e.g.
# quinemc(4222345678921334)) the number of permutations and comparisons grows pretty
# large and can be slow.
Term = namedtuple(
    'Term', 'termset used ones source generation final binary row')

# @coverage


def canonical(item, highorder_a=True, includef=False):
    '''
    Takes an integer and returns the boolean algebra canonical disjunctive normal form.

    Arguments:
    x (integer): number to convert to a canonical boolean expression.

    highorder_a (boolean): When "True" (default) A is the high order bit. When False A is
    the low order bit.

    includef (boolean): When True will append "F(x) = " before the result string

    '''
    if not isinstance(item, int):
        return ValueError(item, "canonical(x) requires an integer.")

    binary = format(item, 'b')
    binary = binary[::-1]
    # No. of letters needed is equal to the length of the binary number representing
    # the length of our number. (E.g. len('11111000') == (8 - 1) == 0b111. len('111') = 3,
    # so we will need A, B, C.
    letters = len(format(len(binary) - 1, 'b'))
    if letters == 1:
        letters = 2

    # Funky formatter: essentially `format('10', '05b')` to get 00010--i.e. a binary equal to
    # the length of letters for each bit position that equals 1 in our input
    indexes = [(format(i.start(), '0' + str(letters) + 'b')) for i in re.finditer('1', binary)]

    miniterms = [_minterms_(m, highorder_a) for m in indexes[::-1]]
    miniterms = sorted(miniterms, reverse=True)

    result = ' + '.join(miniterms)

    if includef is True:
        result = "f(" + str(item) + ") = " + result
    return result

def _minterms_(terms, highorder_a):
    alpha = sorted(string.ascii_letters)
    result = ''
    if highorder_a is False:
        terms = terms[::-1]

    # convert 010 to A'BC'
    for i, terms in enumerate(terms):
        result += alpha[i]
        if terms == '0':
            result += "'"
    return result


# --- END OF Canonical functions ---
# --- START: Go from Minform to Canonical ---
def to_cdnf(min_form, ranged=False):
    """ CONVERT BACK
    --make tuples of current terms [('A', "B'"), ("C'", "D'")] -- terms
    --find greatest letter (D)
    --for each term above create set missing pairs [["C", "C'"], ["D", "D'"]]-->missing_list
    -- create all combinations
        missing_combos = list(itertools.product(*missing_list)
    -- merge them
        final.append([set(term) | set(missing) for missing in missing_combos])
    """
    if isinstance(min_form, str):
        min_form = list(re.split(r"[^a-zA-Z']+", min_form))
    elif isinstance(min_form, list):
        pass
    else:
        return ValueError(min_form, "Not valid input")

    letters = ("".join(min_form)).replace("'", "")
    if ranged:
        last_letter = max(letters)
        first_letter = min(letters)
        letters = [chr(i) for i in range(ord(first_letter), ord(last_letter) + 1)]
    # list of list of letters
    min_form = [re.findall("([A-Za-z]'*)", term_letters) for term_letters in min_form]

    final = []
    for term_letters in min_form:
        missing_letters = set(("".join(term_letters)).replace("'", ""))
        missing_list = list(set(letters) - set(missing_letters))
        missing_list = [[q, q + "'"] for q in missing_list]
        missing_combos = list(itertools.product(*missing_list))
        final.extend(["".join(sorted(set(term_letters) | set(missing)))
                      for missing in missing_combos])

    result = sorted(list(set(final)), reverse=True)
    result = ' + '.join(result)

    return result
# --- END: Go from Minform to Canonical ---

# @coverage


def quinemc(myitem, highorder_a=True, full_results=False):
    '''Short for Quinn-McCluskey algorithm
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    Arguments:
    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms. Can be any
    group of "consistent" letters e.g. "rts + rt's + r't's" is valid.

    highorder_a: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method. May be empty.

    >>> quinemc(248)
    'BC + A'
    >>> quinemc(248, False)
    'C + AB'
    >>> quinemc("ABC' + ABC + AB'C' + AB'C + A'BC")
    'BC + A'
    >>> quinemc(["ABC", "ABC'", "AB'C'", "AB'C", "A'BC"])
    'BC + A'
    >>> result, term_list, possibles = quinemc(248, False, True)

    Invalid input will return a ValueError (e.g. quinemc(["A'BC", "AB"]) is invalid
    because second term must contain a "C".

    Add code for doing dont_care items
    '''
    if (isinstance(myitem, list) and
            len(myitem) == 2 and
            all(isinstance(part, list) for part in myitem)):
        cdnf = _convert_to_terms_(myitem[0], highorder_a)
        if all(isinstance(i, int) for i in myitem[1]):
            dont_care = myitem[1]
        elif all(isinstance(a, str) for a in myitem[1]): # don't care are terms e.g ABC'D
            temp_terms = [set(re.findall("([A-Za-z]'*)", x)) for x in myitem[1]]
            dont_care = [int(_make_binary(z), 2) for z in temp_terms]
    else:
        cdnf = _convert_to_terms_(myitem, highorder_a)
        dont_care = None

    if cdnf is None:
        return ValueError(myitem, "Invalid input")

    test_string = "".join(sorted(re.sub("'", '', cdnf[0])))
    for item in cdnf:
        if test_string != "".join(sorted(re.sub("'", '', item))):
            return ValueError("Term: ", item, " doesn't match valid test ", test_string)

    if full_results:
        return _minimize_(cdnf)
    else:
        return _minimize_(cdnf)[0]

def _convert_to_terms_(item_in, highorder_a):
    result = item_in
    if isinstance(item_in, int):
        result = canonical(item_in, highorder_a).split(' + ')
    elif isinstance(item_in, str):
        result = re.split(r"[^a-zA-Z']+", item_in)
    elif isinstance(item_in, list) and all(isinstance(x, str) for x in item_in):
        result = item_in
    else:
        result = None

    return result

# @coverage
def _minimize_(cdnf):
    """
    Basic driver for performing the over all minimization. Most of the "heavy" work
    is done in _implicants_()

    """
    done = False
    current_generation = 1
    # Step 1: take terms from the canonical form and putting them into generation one
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
    result = sorted(result, reverse=True)
    result = " + ".join(result)

    # Handle special cases-- quinemc(0), quinemc(15), quinemc(255),
    # quinemc("AB + A'B + AB' + A'B'"), etc.
    if term_list and not term_list[0].termset:
        result = "0"
    if result == "":
        result = "1"

    return result, term_list, possibles


# @coverage
def _create_first_generation_(terms):
    my_letters = set(sorted(string.ascii_letters))
    temp_terms = [set(re.findall("([A-Za-z]'*)", x))
                  for x in terms]  # Convert to list of sets

    # Remove duplicate terms if called with something like quinemc("ABCD + CDBA + ABC'D + DC'AB")
    temp_terms = list(temp_terms for temp_terms, _ in itertools.groupby(temp_terms))

    temp_list = [Term(x, False, len(x.intersection(my_letters)), None, 1, None,
                      _make_binary(x), None) for x in temp_terms]
    temp_list = sorted(temp_list, key=attrgetter('ones'))
    for idx, item in enumerate(temp_list):
        temp_list[idx] = item._replace(source=[idx], row=idx)
    return temp_list


def _make_binary(new_term):
    new_term = "".join(sorted(list(new_term)))
    new_term = re.sub("[A-Za-z]'", "0", new_term)
    new_term = re.sub("[A-Za-z]", "1", new_term)
    return new_term


# @coverage
def _merge_terms_(term_list, gen):
    done = False
    new_terms = _create_new_terms_(term_list, gen)

    if new_terms:
        term_list.extend(new_terms)
    else:
        done = True

    return done


# @coverage
# @profile
def _create_new_terms_(orig_term_list, gen):
    # Takes a generation and puts it into a list of lists of terms grouped by the 
    # number of "ones". 
    # Then successively compares items in one list with items from the next list to
    # find terms that can be merged
    used_dict = {}  # a dictionary for used items
    sources = []  # avoid duplicate merges
    result = []

    working_list = [xterms for xterms in orig_term_list if xterms.generation == gen]
    working_list = [list(group) for key, group in
                    itertools.groupby(working_list, attrgetter('ones'))]
    # last_ones = len(working_list) - 1
    current = 0

    while current < len(working_list) - 1:
        # first_list, second_list = working_list[current], working_list[current + 1]
        for xterms in working_list[current]:
            for yterms in working_list[current + 1]:
                sym_set = yterms.termset.symmetric_difference(xterms.termset)
                # doing list(sym_set)[0][0] == list(sym_set)[1][0]: is expensive
                if (len(sym_set) == 2 and
                        sym_set.pop().replace("'", "") == sym_set.pop().replace("'", "")):
                    used_dict[xterms.row] = True
                    used_dict[yterms.row] = True
                    new_term = yterms.termset.intersection(xterms.termset)
                    source = sorted(yterms.source + xterms.source)
                    if source not in sources:
                        sources.append(source)
                        result.append(Term(new_term, False, len(new_term.intersection(
                            set(string.ascii_letters))), source, (gen + 1), None, None, None))
        current += 1
    result = sorted(result, key=attrgetter('ones'))
    for idx, _ in enumerate(result):
        result[idx] = result[idx]._replace(row=len(orig_term_list) + idx)

    # set used terms as used in orig_term_list
    for key in used_dict:
        current = orig_term_list[key]
        orig_term_list[key] = Term(current.termset, True, current.ones,
                                   current.source, current.generation, None,
                                   current.binary, current.row)

    return result


# ---- Finding Implicants / Final Result ---- #
# _implicants_, _required_sources_, _make_find_dict_, and _check_combinations_ do the heavy
# lifting to find the final reduced form
# @coverage
def _implicants_(term_list):
    '''
    Finds terms that will cover unused cases if the essential prime implicants are not
    enough.
    '''
    possible_terms = defaultdict(list)
    list_of_sources = []

    for sources in [zed for zed in term_list if zed.used is False]:
        list_of_sources += list(itertools.chain(sources.source))

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


# @coverage
def _get_columns_(term_list, required):
    """
    term_list -- full list of terms
    required -- integers for terms that are essential prime implicants . . .
        each required int will appear in the source list for only 1 item in needed
    """
    ignore = []
    keep = []

    for index, term in [(i, v) for i, v in enumerate(term_list) if v.used is False]:
        # Find Terms in "needed" that exist in required, add them to the final result,
        # and add that Term's sources to the "columns" we can now ignore (already covered
        # terms)
        if len((set(required) & set(term.source))) >= 1:
            term_list[index] = term._replace(final="Required")
            ignore += itertools.chain(term.source)
        # Otherwise add the sources to our list of "columns" we need to keep
        else:
            keep += itertools.chain(term.source)

    # create a list of the remaining 1st gen terms that we still need to find minterms for
    keep = list(set(keep) - set(ignore))

    return keep


# @coverage
def _make_find_dict_(term_list, keep_columns):
    # Creates a dictionary referencing the remaining tuples that can potentially complete
    # the minimized form.
    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = {}
    for idx, item in [(i, k) for i, k in enumerate(term_list)
                      if k.used is False and k.final is None]:
        temp_source = set(keep_columns) & set(item.source)

        if temp_source:
            temp_tuple = search_tuple(temp_source, len(item.termset))
            find_dict[idx] = temp_tuple

    return find_dict


# @coverage
# @profile
def _check_combinations_(find_dict, term_list, keep_columns):
    matches = []
    possible_terms = defaultdict(list)
    min_length = 0
    break_count = 0

    for fixme in range(2, (len(find_dict) + 1)):
        # adding more and more combinations isnt likely to improve (shorten) length of result
        # so once matches are found we limit how many more sets of combinations we check
        # there may however be funky corner cases
        if break_count >= 2:
            break
        else:
            if matches:
                break_count += 1
        for items in itertools.combinations(find_dict.keys(), fixme):
            combined_sources = set()
            temp_count = 0
            for idx in items:
                combined_sources.update(find_dict[idx].sourceSet)
                temp_count += find_dict[idx].length
            if (set(keep_columns) == combined_sources
                    and (min_length == 0 or temp_count <= min_length)):
                if (temp_count < min_length or min_length == 0):
                    del matches[:]
                    min_length = temp_count
                matches.append(items)

    if matches:
        for idx, value in enumerate(matches):
            for i in value:
                if idx == 0:
                    term_list[i] = term_list[i]._replace(final="Added")
                possible_terms[idx].append(term_list[i])

    return possible_terms


if __name__ == "__main__":
    #quinemc(42589768824798729982179, True, True)
    #print(canonical(9927465))
    #quinemc(638, 1, 1)
    quinemc(to_cdnf("A + C", 1))
    #to_cdnf("B'CD + A'C'D' + A'B'D'")
    #quinemc(canonical(27856))
    # quinemc(to_cdnf("A + FI"))
    #import doctest
    # doctest.testmod()
    #print("2003", quinemc(2003))
    #print("2077", quinemc(2077))
    #print("2078", quinemc(2078))
    #print("2046", quinemc(2046))
    #print("255", quinemc(255))
    #print("0", quinemc(0))
    # canonical(2046)
    pass

# 65024
# 15872
# 48640
# **Situation 4 -- 38400**
# interesting one: a, b, c = quinemc(248725692, True, True)
