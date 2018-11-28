""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
import utils
from numpy import zeros, ndarray
from shared import *
from evaluation import index_isoform_locations

ALPHABET = [TERMINATOR] + BASES

def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
    
    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """
    suffixes = utils.get_suffixes(s)
    suffixes.sort()
    return [suf[1] for suf in suffixes if suf is not None]

def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    L = ''
    for i in sa:
        L += s[i-1]
    return L

def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted
    """
    return ''.join(sorted(L))

def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    """
    M = {}
    for ch in ALPHABET:
        try:
            M[ch] = F.index(ch)
        except ValueError:
            M[ch] = -1
    return M

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    OCC = {ch: ([0]*len(L)) for ch in ALPHABET}
    
    for index in range(len(L)):
        for let in range(len(ALPHABET)):

            if L[index] == ALPHABET[let]:
                is_same = 1
            else:
                is_same = 0    

            OCC[ALPHABET[let]][index] = OCC[ALPHABET[let]][index-1] + is_same

    return OCC        

def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    """
    _length = 0
    _range = None
    i = len(p) - 1

    sp = M[p[i]]
    if sp == -1:
        return (_range, _length)

    ep = sp + occ[p[i]][-1] - 1
    i -= 1
    _length += 1
    _range = (sp, ep)

    while i >= 0 and sp <= ep:

        x = (occ[p[i]][sp - 1])
        y = (occ[p[i]][ep])

        sp = (M[p[i]] + x)
        ep = (M[p[i]] + y) - 1
        
        if sp > ep:
            break
        else:    
            i -= 1
            _length += 1
            _range = (sp, ep)

    _range = (_range[0], _range[1] + 1)


    return (_range, _length)    

      

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        
        self.transcriptome, self.index_dict = self.build_transcriptome(known_genes, genome_sequence)
        #print('Built transcriptome')
        sequence_to_search = self.transcriptome  # or genome_sequence?

        self._sa = get_suffix_array(sequence_to_search + '$')
        #print('Built SA')
        self._bwt = get_bwt(sequence_to_search + '$', self._sa)
        #print('Built BWT')
        self._F = get_F(self._bwt)
        self._M = get_M(self._F)
        self._occ = get_occ(self._bwt)
        #print('Built OCC matrix')

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        match_sequence = read_sequence
        bases = {'A', 'C', 'G', 'T'}
        best_len = 0
        curr_length = 0

        range_and_len = exact_suffix_matches(match_sequence, self._M, self._occ)
        rng, length = range_and_len
        match_rng = rng
        curr_length = length + 1

        while curr_length < len(read_sequence):
            
            best_len = 0
            swap_bases = bases.difference(match_sequence[-curr_length])
            #print('$$$$')

            for ch in swap_bases:
                new_match = match_sequence[:-curr_length] + ch + match_sequence[-curr_length+1:]
                #print(new_match)
                rng, length = exact_suffix_matches(new_match, self._M, self._occ)
                if length > best_len:
                    match_sequence = new_match
                    best_len = length
                    match_rng = rng

            curr_length = best_len + 1


        return self.get_match_locations(read_sequence, match_sequence, self._sa[match_rng[0]])  #index is location in the transcriptome
       


    def get_match_locations(self, original, actual, index):
        ctr = 0
        previous = -2
        match_locations = []
        curr_match = (None, None, None)  #(read, ref, len)


        for a,b in zip(original, actual):
            #print(a, b)
            if a == b:
                ref_index = self.index_dict[index+ctr]
                if curr_match[0] == None:
                    curr_match = (ctr, ref_index, 1)
                    #print(curr_match)
                elif ref_index - previous != 1:
                    #print('Jumped intron')
                    match_locations.append(curr_match)
                    curr_match = (ctr, ref_index, 1)
                else:
                    #print('Adding 1')
                    curr_match = (curr_match[0], curr_match[1], curr_match[2] + 1)

                previous = ref_index    

            else:
                match_locations.append(curr_match)
                curr_match = (None, None, None)

                previous = -2

            ctr += 1

        if curr_match[0] != None:
            match_locations.append(curr_match)

        return match_locations    
             


    def build_transcriptome(self, known_genes, genome_sequence):
        transcriptome = ''
        index_dict = {}
        length = 0
        isoforms = {}
        seq = ''
        for g in known_genes:
            for i in g.isoforms:
                for e in i.exons:
                    segment = genome_sequence[e.start:e.end]
                    seq += segment
                    for x in range(e.end - e.start):
                        index_dict[length + x] = e.start + x
                    transcriptome += segment
                    length += len(segment)
                isoforms[i.id] = seq
                seq = ''
                transcriptome += '!' #Assuming reads can not span across isoforms so they are seperated with unique char
                length += 1
        return transcriptome, index_dict