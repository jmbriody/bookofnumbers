import re
import math
import itertools
from collections import namedtuple, defaultdict
from operator import attrgetter
from profilehooks import coverage, profile

Term = namedtuple('Term', 'termset used ones source generation')

def canonical(x, INCLUDEF=True):
    """
    Takes an integer and returns the boolean algebra canonical disjunctive normal form.

    """
    binary = format(x, 'b') # obviously convert to binary
    letters = _letters_(len(binary)) # determine the number of letters a minterm will contain
    # Reverse the binary and create a list of all the positions with a "1".
    binary = binary[::-1]
    indexes = [(i.start() + 1) for i in re.finditer('1', binary)]
    # Use the position list to generate our minterms
    miniterms = [_miniterms_(m, letters) for m in indexes[::-1]]

    result = ' + '.join(miniterms)
    if INCLUDEF is True:
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
    would have 3 letters representing 256 possible combinations.
    """
    letters = 2
    start = 4
    while length > start:
        letters = letters + 1
        start = start * 2

    return letters


def _miniterms_(m, letters, REVERSE=False):
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

https://www.allaboutcircuits.com/textbook/digital/chpt-8/logic-simplification-karnaugh-maps/

http://www.plantation-productions.com/Webster/www.artofasm.com/Linux/PDFs/DigitalDesign.pdf
    (Section 3.4)

    """
    alpha = 'ABCDEFG'
    current_letter = 0
    control = 1
    result = ''
    while letters > 0:
        # first add the letter because we know we need it. If REVERSE is true the last
        # letter will be 01010101, 2nd to last will be 00110011, etc.
        if REVERSE is True:
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

def qmc(nIn):
    result = quinmc(nIn)

    return result[0]

# @coverage
def quinmc(nIn):
    """
    https://www.allaboutcircuits.com/technical-articles/everything-about-the-quine-mccluskey-method/
    https://pythontips.com/2015/06/06/why-should-you-use-namedtuple-instead-of-a-tuple/
    z = quinmc(4222345678921334)
    """
    cdnf = canonical(nIn, False).split(' + ')
    term_list = _convert_to_tuples_(cdnf, 1)
    DONE = False
    current_generation = 1

    while not DONE:
        term_list, DONE = _merge_terms_(term_list, current_generation)
        current_generation += 1

    final_result, possibles = _minimize_(nIn, [x for x in term_list if x.used is False])

    result = ["".join(sorted(tempItem.termset)) for tempItem in final_result]

    result = " + ".join(result)

    if nIn == 0:
        result = "0"
    elif result == "":
        result = "1"

    return result, final_result, term_list, possibles

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

@coverage
def _minimize_(cdnum, needed):
    # find essential prime implicants
    required = _required_sources_(needed)
    needed, final_list, keep_columns = _get_columns_(needed, required)

    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = {}

    for idx, item in enumerate(needed):
        temp_source = set(keep_columns) & set(item.source)

        if len(temp_source) > 0:
            temp_tuple = search_tuple(temp_source, len(item.termset))
            find_dict[idx] = temp_tuple

    combos = len(find_dict) + 1
    FINISHED = False
    key_list = [x for x in find_dict.keys()]
    possible_terms = defaultdict(list)
    matches = []

    while not FINISHED:
        # check if single term will "cover" remaining items e.g. qmc(2077)
        for idx, val in find_dict.items():
            if set(keep_columns) == val.sourceSet:
                print(cdnum)
                final_list.append(needed[idx])
                FINISHED = True
                break
        if FINISHED:
            break

        for x in range(2, combos):
            # currentCombos = list(itertools.combinations(key_list, x))
            for items in itertools.combinations(key_list, x):
                combined_sources = set()
                for idx in items:
                    combined_sources.update(find_dict[idx].sourceSet)
                if set(keep_columns) == combined_sources:
                    matches.append(items)
            if len(matches) > 0:
                FINISHED = True
                for idx, value in enumerate(matches):
                    for i in value:
                        if idx == 0:
                            final_list.append(needed[i])
                        possible_terms[idx].append(needed[i])
                break
        FINISHED = True

    return final_list, possible_terms


def _convert_to_tuples_(terms, gen, source=None):
    my_letters = set('ABCDEFG')
    temp_terms = [set(re.findall("([A-Z]'*)", x)) for x in terms] # Convert to list of sets
    temp_list = [Term(x, False, len(x.intersection(my_letters)), source, gen) for x in temp_terms]
    temp_list = sorted(temp_list, key=attrgetter('ones'))
    for idx, item in enumerate(temp_list):
        temp_list[idx] = item._replace(source=[idx])
    return temp_list


def _merge_terms_(orig_term_list, gen):
    working_list = [x for x in orig_term_list if x.generation == gen]
    used_dict = {} # a dictionary for used items
    new_tuples = []
    sources = []    # avoid duplicate merges
    DONE = False

    for x in working_list:
        for y in working_list:
            if y.ones == (x.ones + 1):
                # foobared my logic--band aid fix
                sym_set = y.termset.symmetric_difference(x.termset)

                if len(sym_set) == 2 and list(sym_set)[0][0] == list(sym_set)[1][0]:
                    x_index, y_index = orig_term_list.index(x), orig_term_list.index(y)
                    used_dict[x_index] = True
                    used_dict[y_index] = True
                    new_term = y.termset.intersection(x.termset)
                    source = sorted(y.source + x.source)
                    if source not in sources:
                        sources.append(source)
                        new_tuples.append(Term(new_term, False, len(new_term.intersection('ABCDEFG')), source, (gen + 1)))
    new_tuples = sorted(new_tuples, key=attrgetter('ones'))
    for key in used_dict.keys():
        current = orig_term_list[key]
        orig_term_list[key] = Term(current.termset, True, current.ones, current.source, current.generation)
    if len(new_tuples) > 0:
        orig_term_list.extend(new_tuples)
    else:
        DONE = True

    return orig_term_list, DONE

if __name__ == "__main__":
    mj = canonical(14)
    mk = canonical(248)
    mi = canonical(22256)
    print(mj, mk, mi)
    for q in range(2000, 2100):
        print(q, qmc(q))
