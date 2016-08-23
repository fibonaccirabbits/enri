'''
Predicts and writes an output file for top n predicted conformations. The
output file is writen to the input directory

INPUT: desc_merged.txt
OUTPUT: *predicted*.txt
ARGUMENTS: input_path, number of desired output (n), beta, ranker (wp or p)
USAGE: python descriptors2predictions desc_merged.txt, n, beta,ranker

'''

#import stuff
from enri import Enri
import sys

e = Enri()
input_path = sys.argv[1]
n = int(sys.argv[2])
beta = float(sys.argv[3])
ranker = sys.argv[4]
path, name, name_ext = e.path2names(input_path)
outdir = path
e.file2top_predicted2(input_path, n, beta, ranker, outdir)
