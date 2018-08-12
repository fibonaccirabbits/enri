### Read the paper [here](https://onlinelibrary.wiley.com/doi/full/10.1111/cbdd.12900)

ENRI
----
ENRI is a tool for selecting structure-based virtual screening targets. The
tool is a binary classifier coupled with a synthetic over-sampling
procedure. ENRI, currently, comprises four programs: enri.py,
pdb2descriptors.py, descriptors2predictions.py, plot_hist.py.

PREQUISITES
-----------
A fully funtional DogSiteSCorer.
A python plotting library: matplotlib.*
A python tabulation library.**

*only when plotting is desired
**only when tabulation is desired

enri.py
------
This is the main program where all of ENRI's functionalities are defined.

pdb2descriptors.py
------------------
Extracts pockets and descriptors from pdb files. Interfaces with DoGSiteScorer.
Please make sure you have fully funtional DoGSiteScorer.

INPUT: pdb_path
OUTPUT: desc_merged.txt
ARGUMENTS: pdb_path
USAGE: python pdb2descriptors.py pdb_path
EXAMPLE: python pdb2descriptors.py /enri_rc8/sample_files/pdbdir

descriptors2predictions.py
--------------------------

Predicts and writes an output file for top n predicted conformations. The
output file is written to the input directory

INPUT: desc_merged.txt
OUTPUT: *predicted*.txt
ARGUMENTS: input_path, number of desired output (n), over-sampling paramter (beta), ranker (wp or p)
USAGE: python descriptors2predictions desc_merged.txt, n, beta,ranker
EXAMPLE: python descriptors2predictions.py /enri_rc8/sample_files/descdir/desc_merged.txt 10 0.5 wp

plot_hist.py
------------
Plots histogram from a desc_merged.txt file.

INPUT: file_path
OUTPUT: *descriptorname.pdf
ARGUMENTS: file_path
USAGE: python plot_hist.py file_path
EXAMPLE: python plot_hist.py /enri_rc8/sample_files/descdir/desc_merged.txt
 
