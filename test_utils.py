import utils
from project import *

stringg = 'ATCG$'

suf = utils.get_suffixes(stringg)



"""
final = utils.sort_suffixes(suf, 4)
#bomb = [s[1] for s in final if s is not None]
for f in final:
    print(f)
"""    
print(get_occ('ATG$AAAGTG'))    