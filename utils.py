""" 
    Helper functions for building suffix arrays.
"""




def get_suffixes(s):
    """
    Given a string s, returns a list of tuples of the form (suffix, index)
    """
    suffixes = []
    for i in range(len(s)):
        suffixes.append((s[i:] + '$'*(i), i))

    return suffixes

def sort_suffixes(suffixes, depth):
    """
    Sorts the given suffixes (tuple of the form (suffix, index)) using only the first depth characters of each suffix. 
    Returns a sorted list of tuples (suffix, index).
    The list of suffix indices will be sorted with the order ($, A, C, G, T)
    """
    if len(suffixes) == 0:
        return [None]

    if len(suffixes) == 1 or depth <= 0:
        return suffixes

    bins = {c:[] for c in '$ACGT'}
    for s in suffixes:
        bins[s[0][0]] += [(s[0][1:], s[1])]

    #big_one = [sort_suffixes(bins[ch], depth-1) for ch in '$ACGT' if bins[ch] is not None]
    

    return (sort_suffixes(bins['$'], depth-1) 
        + sort_suffixes(bins['A'], depth-1) 
        + sort_suffixes(bins['C'], depth-1) 
        + sort_suffixes(bins['G'], depth-1) 
        + sort_suffixes(bins['T'], depth-1))
    """
    (sort_suffixes(bins['$'], depth-1) 
        + sort_suffixes(bins['A'], depth-1) 
        + sort_suffixes(bins['C'], depth-1) 
        + sort_suffixes(bins['G'], depth-1) 
        + sort_suffixes(bins['T'], depth-1))
    """    