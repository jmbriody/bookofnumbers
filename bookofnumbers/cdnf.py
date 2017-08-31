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

In the standard mode (reverse = False) A is the low order bit:

    BP  |  C  |  B  |  A     248 (11111000)
    ---------------------    ---
    1   |  0  |  0  |  0      0
    2   |  0  |  0  |  1      0
    3   |  0  |  1  |  0      0
    4   |  0  |  1  |  1      1     ABC'
    5   |  1  |  0  |  0      1     A'B'C
    6   |  1  |  0  |  1      1     AB'C
    7   |  1  |  1  |  0      1     A'BC
    8   |  1  |  1  |  1      1     ABC

When reverse is True C would be the low order bit for 248.

Two primary user functions are:
caononical(x, reverse=False, includef=False):
    Simplest call is canonical(248) which will return the canonical boolean algebraic expression
    for "248". In the default setup A is the Low Order bit (reverse=False)--i.e. 0101. When reverse
    is True A will be the High Order bit--i.e. 0000111100001111. includef determines whether
    "F(248) = " should be prepended to the result.



>>> canonical(248)
"ABC + A'BC + AB'C + A'B'C + ABC'"
>>> canonical(248, True)
"CBA + C'BA + CB'A + C'B'A + CBA'"
>>> canonical(248, False, True)
"f(248) = ABC + A'BC + AB'C + A'B'C + ABC'"

quinmc(myitem, reverse=False, full_results=False):
    Short for Quinn-McCluskey algorithm:
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms

    reverse: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --final_result: List of Named Tuples "termset" containing terms of the Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method. May be empty.

>>> quinmc(248)
'AB + C'
>>> quinmc(248, True)
'BC + A'
>>> quinmc("ABC + A'BC + AB'C + A'B'C + ABC'")
'AB + C'
>>> quinmc(["ABC", "A'BC", "AB'C", "A'B'C", "ABC'"])
'AB + C'
>>> result, final_result, term_list, possibles = quinmc(248, False, True)

"""

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
# quinmc(4222345678921334)) the number of permutations and comparisons grows pretty
# large and can be slow.
Term = namedtuple('Term', 'termset used ones source generation')

def canonical(x, reverse=False, includef=False):
    """
    Takes an integer and returns the boolean algebra canonical disjunctive normal form.

    Arguments:
    x: integer to convert to a canonical boolean expression. based on Section 3.4 of
    "Digital Design" from Randall Hyde's "The Art of Assembly Language Programming"
http://www.plantation-productions.com/Webster/www.artofasm.com/Linux/PDFs/DigitalDesign.pdf

    reverse: When "False" (default) A is the low order bit. When True A is the high order
    bit.

    includef: When True will append "F(x) = " before the result string

    """
    binary = format(x, 'b') # obviously convert to binary
    letters = _letters_(len(binary)) # determine the number of letters a minterm will contain
    # Reverse the binary and create a list of all the positions with a "1".
    binary = binary[::-1]
    indexes = [(i.start() + 1) for i in re.finditer('1', binary)]
    # Use the position list to generate our minterms
    miniterms = [_miniterms_(m, letters, reverse) for m in indexes[::-1]]

    result = ' + '.join(miniterms)
    if includef is True:
        result = "f(" + str(x) + ") = " + result
    return result


def _letters_(length):
    """
    Based on the length of the binary number determine the number of letters
    in a minterm.

    We start with 2 letters ( A and B ) for four bits. Available minterms would be
    A'B', AB', A'B, AB. They can be combined in 16 different ways when taking any
    combination with 1 to 4 minterms. Essentially four bits can represent 16 different
    equations.

    We add 1 letter each time we double the number of bits. Anything from 5 to 8 bits
    would have 3 letters representing 256 possible combinations etc.
    """
    letters = 2
    start = 4
    while length > start:
        letters = letters + 1
        start = start * 2

    return letters


def _miniterms_(m, letters, reverse=False):
    """
    Generate a minterm based on the location/index of the binary digit (m) and the number of letters

    Reference table:

    BP -- Bit Position

    BP  |  C  |  B  |  A
    ---------------------
    1   |  0  |  0  |  0
    2   |  0  |  0  |  1   <-- ceiling(BP/1) % 2 != 0 [A]
    3   |  0  |  1  |  0       <-- ceiling(BP/2) % 2 != 0 [B]
    4   |  0  |  1  |  1
    5   |  1  |  0  |  0           <-- ceiling(BP/4) % 2 != 0 [C]
    6   |  1  |  0  |  1
    7   |  1  |  1  |  0
    8   |  1  |  1  |  1
                                        <-- ceiling(BP/8) % 2 != 0 [D]

    For example BP == 7.
    [A] ceiling(7/1) % 2 != 0 --> True -- So we add "'" to A
    [B] ceiling(7/2) % 2 != 0
        4 % 2 != 0 --> False -- B is a "1" so it is unprimed
    [C] ceiling(7/4) % 2 != 0
        2 % 2 != 0 --> False -- C is a "1" so it is unprimed

    So we get A'BC for BP 7 or C'BA if REVERSE it True

    """
    alpha = sorted(string.ascii_letters)
    current_letter = 0
    control = 1
    result = ''

    while letters > 0:
        # first add the letter because we know we need it. If REVERSE is true the last
        # letter will be 01010101, 2nd to last will be 00110011, etc.
        if reverse is True:
            result = result + alpha[letters - 1]
        else:
            result = result + alpha[current_letter]
        # See reference table.
        if math.ceil(m / control) % 2 != 0:
            result = result + "'"
        control = control * 2
        current_letter = current_letter + 1
        letters = letters - 1
    return result
# --- END OF Canonical functions ---

# @coverage
def quinmc(myitem, reverse=False, full_results=False):
    '''Short for Quinn-McCluskey algorithm
    https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm

    Arguments:
    myitem: can be an integer, string representing a canonical boolean expression, or a list
    of minterms for a canonical expression.
    --integer: first passed to canonical to get the normal form
    --string or list: must contain minterms that are canonical normal form terms

    reverse: Same as for canonical. Only has any affect if "myitem" is an int.

    full_results: In the default mode "False" only a string containing the minimized function is
    returned. When True a list of various data structures is returned.
    --result: Minimized expression
    --final_result: List of Named Tuples "termset" containing terms of the Minimized expression
    --term_list: List of Named Tuples "termset" of all terms used in the reduction/minimize process
    --possibles: A dictionary of possible minterms when the essential prime implicants do not cover
    the whole canonical form. Essentially items generated using Petrick's Method. May be empty.

    >>> quinmc(248)
    'AB + C'
    >>> quinmc(248, True)
    'BC + A'
    >>> quinmc("ABC + A'BC + AB'C + A'B'C + ABC'")
    'AB + C'
    >>> quinmc(["ABC", "A'BC", "AB'C", "A'B'C", "ABC'"])
    'AB + C'
    >>> result, final_result, term_list, possibles = quinmc(248, False, True)

    Invalid input will return a ValueError (e.g. quinmc(["A'BC", "AB"]) is invalid
    because second term must contain a "C".
    '''
    if isinstance(myitem, int):
        cdnf = canonical(myitem, reverse).split(' + ')
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
        return minimize(cdnf)
    else:
        return minimize(cdnf)[0]

# @coverage
def minimize(cdnf):
    """
    Basic driver for performing the over all minimization. Most of the "heavy" work
    is done in _implicants_()

    """
    # Step 1: take terms from the caononical form and putting them into generation one
    # of our term_list
    term_list = _convert_to_tuples_(cdnf, 1)
    done = False
    current_generation = 1

    # Step 2: merge terms of each generation to create next generation until no more merges
    # are possible (_merge_terms_ and _create_new_tuples_)
    while not done:
        term_list, done = _merge_terms_(term_list, current_generation)
        current_generation += 1

    # Step 3: Generate our final result from all terms in term_list that have not been used in
    # a merge
    final_result, possibles = _implicants_([x for x in term_list if x.used is False])

    result = ["".join(sorted(tempItem.termset)) for tempItem in final_result]
    result = " + ".join(result)

    # Handle special cases-- quinmc(0), quinmc(15), quinmc(255),
    # quinmc("AB + A'B + AB' + A'B'"), etc.
    if len(term_list) == 1 and len(term_list[0].termset) == 0:
        result = "0"
    if result == "":
        result = "1"

    return result, final_result, term_list, possibles

# @coverage
def _convert_to_tuples_(terms, gen, source=None):
    my_letters = set(sorted(string.ascii_letters))
    temp_terms = [set(re.findall("([A-Za-z]'*)", x)) for x in terms] # Convert to list of sets
    temp_list = [Term(x, False, len(x.intersection(my_letters)), source, gen) for x in temp_terms]
    temp_list = sorted(temp_list, key=attrgetter('ones'))
    for idx, item in enumerate(temp_list):
        temp_list[idx] = item._replace(source=[idx])
    return temp_list

# @coverage
def _merge_terms_(orig_term_list, gen):
    working_list = [x for x in orig_term_list if x.generation == gen]
    done = False
    new_tuples, used_dict = _create_new_tuples_(working_list, orig_term_list, gen)

    for key in used_dict.keys():
        current = orig_term_list[key]
        orig_term_list[key] = Term(current.termset, True, current.ones,
                                   current.source, current.generation)
    if len(new_tuples) > 0:
        orig_term_list.extend(new_tuples)
    else:
        done = True

    return orig_term_list, done

# @coverage
def _create_new_tuples_(working_list, orig_term_list, gen):
    used_dict = {} # a dictionary for used items
    sources = []    # avoid duplicate merges
    result = []
    my_letters = set(string.ascii_letters)
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
                        result.append(Term(new_term, False, len(new_term.intersection(my_letters)),
                                           source, (gen + 1)))
    result = sorted(result, key=attrgetter('ones'))

    return result, used_dict
    
# ---- Finding Implicants / Final Result ---- #
# _implicants_, _required_sources_, _make_find_dict_, and _check_combinations_ do the heavy
# lifting to find the final reduced form
# @coverage
def _implicants_(needed):
    '''
    Finds terms that will covered unused cases if the essential prime implicants are not
    enough.
    '''
    possible_terms = defaultdict(list)
    required = _required_sources_(needed)
    needed, final_list, keep_columns = _get_columns_(needed, required)

    # if _get_columns_ ends with nothing in keep_columns it means essential prime implicants
    # are all that is needed so we are done
    if len(keep_columns) == 0:
        finished = True
    else:
        finished = False

    # check if single term will "cover" remaining items e.g. qmc(2077)
    if not finished:
        find_dict = _make_find_dict_(needed, keep_columns)
        for idx, val in find_dict.items():
            if set(keep_columns) == val.sourceSet:
                final_list.append(needed[idx])
                finished = True
                break

    # if a single term doesn't cover the remaining 1st gen items start looking
    # for combinations of Terms that will fit the bill.
    if not finished:
        final_list, possible_terms = _check_combinations_(find_dict,
                                                          final_list, keep_columns, needed)

    return final_list, possible_terms

# @coverage
def _required_sources_(needed):
    """
    input: needed -- a list of "Term" namedtuples where used==False
    output: a list of integers for "Term" items from generation 1 that
        appear only once

    In a nutshell this finds reduced terms that are essential prime implicants--items
    that have an ancestor that has been used only once.
    """
    list_of_sources = []

    for x in needed:
        list_of_sources += list(itertools.chain(x.source))

    required = [x for x in list_of_sources if list_of_sources.count(x) == 1]

    return required

# @coverage
def _get_columns_(needed, required):
    """
    needed -- list of "Term" tuples that have used==False
    required -- integers for terms that are essential prime implicants . . .
        each required int will appear in the source list for only 1 item in needed
    """
    final_l = []
    ignore = []
    keep = []

    for x in needed:
        # Find Terms in "needed" that exist in required, add them to the final result,
        # and add that Term's sources to the "columns" we can now ignore (already covered
        # terms)
        if len((set(required) & set(x.source))) == 1:
            final_l.append(x)
            ignore += set(itertools.chain(x.source))
        # Otherwise add the sources to our list of "columns" we need to keep
        else:
            keep += set(itertools.chain(x.source))

    # create a list of the remaining 1st gen terms that we still need to find minterms for
    keep = list(set(keep) - set(ignore))

    # Clean up the "needed" list by removing the terms we've used as essential prime implicants
    for x in final_l:
        needed.remove(x)

    return needed, final_l, keep

# @coverage
def _make_find_dict_(needed, keep_columns):
    # Creates a dictionary referencing the remaining tuples that can potentially complete
    # the minimized form. 
    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = {}
    for idx, item in enumerate(needed):
        temp_source = set(keep_columns) & set(item.source)

        if len(temp_source) > 0:
            temp_tuple = search_tuple(temp_source, len(item.termset))
            find_dict[idx] = temp_tuple
    return find_dict

# @coverage
def _check_combinations_(find_dict, final_list, keep_columns, needed):
    # key_list = [x for x in find_dict.keys()]
    match_idx = []
    matches = []
    possible_terms = defaultdict(list)
    done = False

    for x in range(2, (len(find_dict) + 1)):
        for items in itertools.combinations(find_dict.keys(), x):
            if done:
                break
            combined_sources = set()
            temp_count = 0
            for idx in items:
                combined_sources.update(find_dict[idx].sourceSet)
                temp_count += find_dict[idx].length
            if set(keep_columns) == combined_sources:
                matches.append(items)
                match_idx.append(temp_count)        # find result with fewest literals
        # Found at least 1 set of terms that will work
        if len(matches) > 0:
            min_index = match_idx.index(min(match_idx))
            for idx, value in enumerate(matches):
                for i in value:
                    if idx == min_index:
                        final_list.append(needed[i])
                    # also create a list of all options that will fit the bill
                    possible_terms[idx].append(needed[i])
                    done = True
            break
    return final_list, possible_terms





if __name__ == "__main__":
#    print("here")
#    a, b, c, d = quinmc(42589768824798729982179, False, True)
#    print(a)
#    a = quinmc(2077)
    for x in range(2000,  2100):
        print(x)
        quinmc(x)
        
#    a,  b,  c,  d = quinmc(2003, False, True)
#    for i,  v in d.items():
#        print(i, v)
#    mj = canonical(14)
# skokie--6:30 -- 8:45 concert
#    mk = canonical(248)
#    mi = canonical(22256)
#    print(mj, mk, mi)

#    a, b, c, d = quinmc(2078)
#    print(a)
#    for x in b:
#        print(x)
#    print("Termlist")
#    for x in c:
#        print(x)
#    print("Dictionary")
#    for k, v, in d.items():
#        print(k, ":\t", v)
