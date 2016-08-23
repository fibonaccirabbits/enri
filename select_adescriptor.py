'''
writes a descriptor file based on a certain selection kewords arguments
output is written to outfiles directory

INPUT: desc_merged.txt
OUTPUT: *variable*.txt
ARGUMENTS: input_path, pocket_name, variable_name
USAGE: python select_descriptors.py path/to/train3.txt P_1 volume

'''

#import stuff
from enri import Enri
import sys

e = Enri()
input_path = sys.argv[1]
pocket = sys.argv[2]
variable = sys.argv[3]
e.select_adescriptor(input_path, pocket, variable)
