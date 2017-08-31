import string
from collections import namedtuple, defaultdict
from operator import attrgetter
import cdnf
from cdnf import quinmc, _required_sources_, _get_columns_

Term = namedtuple('Term', 'termset used ones source generation')

def _make_find_dict_(needed, keep_columns):
    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = {}
    for idx, item in enumerate(needed):
        temp_source = set(keep_columns) & set(item.source)

        if len(temp_source) > 0:
            temp_tuple = search_tuple(temp_source, len(item.termset))
            find_dict[idx] = temp_tuple
    return find_dict

def _check_combinations_(combos, find_dict, final_list):
    key_list = [x for x in find_dict.keys()]
    match_idx = []
    matches = []
    possible_terms = defaultdict(list)

    for x in range(2, combos):
        for items in itertools.combinations(key_list, x):
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
            finished = True
            for idx, value in enumerate(matches):
                for i in value:
                    if idx == min_index:
                        final_list.append(needed[i])
                    # also create a list of all options that will fit the bill
                    possible_terms[idx].append(needed[i])
            break
    return final_list, possible_terms


def _new_imp_(needed):
    # find essential prime implicants
    required = _required_sources_(needed)
    needed, final_list, keep_columns = _get_columns_(needed, required)

    search_tuple = namedtuple('search_tuple', 'sourceSet length')
    find_dict = _make_find_dict_(needed, keep_columns)
    possible_terms = defaultdict(list)

    finished = False

    # check if single term will "cover" remaining items e.g. qmc(2077)
    for idx, val in find_dict.items():
        if set(keep_columns) == val.sourceSet:
            final_list.append(needed[idx])
            finished = True
            break

    if not finished:
        final_list, possible_terms = _check_combinations_((len(find_dict) + 1), find_dict, final_list)
    
    # if a single term doesn't cover the remaining 1st gen items start looking
    # for combinations of Terms that will fit the bill.
    
    return final_list, possible_terms

zed = cdnf
cdnf._implicants_ = _new_imp_

