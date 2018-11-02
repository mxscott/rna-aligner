import utils
from project import *

string = 'CATCATA$'
p = 'CAT'

M = get_M('$AAACCTT')

occ = get_occ('ATCCT$AA')


print(exact_suffix_matches(p, M, occ))
