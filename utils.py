""" 
    Helper functions for building suffix arrays.
"""
import shared



def get_suffixes(s):
    """
    Given a string s, returns a list of tuples of the form (suffix, index)
    """
    suffixes = []
    for i in range(len(s)):
        suffixes.append((s[i:(i+25)] + '$'*(i), i))

    return suffixes

def sort_suffixes(suffixes, depth):
    """
    NOT USING THIS ANYMORE

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



def parse_tab_file(filepath):
    genes = set()
    curr_gene = ''
    isoforms = []
    curr_iso = ''
    exons = []
    with open(filepath) as f:
        line = f.readline()
        while line != '':
            info = line.strip().split('\t')

            if 'gene' in info[0]:
                if len(exons) != 0:
                    curr_iso = shared.Isoform(curr_iso, exons)
                    isoforms.append(curr_iso)
                    exons.clear()
                if len(isoforms) != 0:
                    curr_gene = shared.Gene(curr_gene, isoforms)
                    genes.add(curr_gene)
                    isoforms.clear()

                curr_gene = info[1]

            elif 'isoform' in info[0]:
                curr_iso = info[1]
                if len(exons) != 0:
                    curr_iso = shared.Isoform(curr_iso, exons)
                    isoforms.append(curr_iso)
                    exons.clear()
                    curr_iso = info[1]
            else:        
                exon = shared.Exon(info[1], int(info[2]), int(info[3]))
                exons.append(exon)
                
            line = f.readline()

    curr_iso = shared.Isoform(curr_iso, exons)
    isoforms.append(curr_iso)
    curr_gene = shared.Gene(curr_gene, isoforms)
    genes.add(curr_gene)

    return genes