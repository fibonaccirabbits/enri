'''
Extracts pockets and descriptors from pdb files. Interfaces with DoGSiteScorer.
Please make sure you have fully funtional DoGSiteScorer.

INPUT: pdb_path
OUTPUT: desc_merged.txt
ARGUMENTS: pdb_path
USAGE: python pdb2descriptors.py pdb_path

'''

#import stuff
import sys
from enri import Enri

e = Enri()
e.pdb2desc_from_path(sys.argv[1])
e.name2firstcol_from_path(sys.argv[1])
e.merge_edt_from_path(sys.argv[1])
