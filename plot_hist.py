'''
Plots histogram from a desc_merged.txt file.

INPUT: file_path
OUTPUT: *descriptorname.pdf
ARGUMENTS: file_path
USAGE: python pdb2descriptors.py file_path

'''

#import stuff
from enri import Enri
import sys


e = Enri()
filepath = sys.argv[1]
path, name, name_ext = e.path2names(filepath)
data, headers = e.parse_desc_merged_txt(filepath)
p1, p1headers = e.get_pocket2(data,headers)
e.plot_hist_nolabel(p1, p1headers,path)
