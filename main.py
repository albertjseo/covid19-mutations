from Bio import SeqIO
import requests
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# load the sequence files
gb_ref_file = "/Users/albertseo/Development/covid19-mutations/ncbi_dataset/data/GenBank/NC_045512.2.gb"
gb_var_file = "/Users/albertseo/Development/covid19-mutations/ncbi_dataset/data/GenBank/MZ571142.1.gb"

# parse the files
ref_seq = SeqIO.read(gb_ref_file, "genbank").seq
var_seq = SeqIO.read(gb_var_file, "genbank").seq

# checkpoint
print(ref_seq)
print(var_seq)

# perform alignment using aligner object
# PairwiseAligner provides more flexibility and better control over alignment parameters
aligner = PairwiseAligner() # establish the object
aligner.mode = "global"  # global (Needleman-Wunsch) alignment
aligner.match_score = 1 # a score that indicates how well a read nucleotide matches a reference nucleotide
aligner.mismatch_score = -1 # the negative penalty assigned when two characters from different sequences are aligned that are not identical
aligner.open_gap_score = -2 # the penalty value assigned to the act of initiating a gap (inserting a space) when aligning two sequences
aligner.extend_gap_score = -0.5 # the penalty assigned to extending a gap by one base in a sequence alignment

alignments = aligner.align(str(ref_seq), str(var_seq))
alignment = alignments[0] # access the best alignment

print("Best Alignment:")
print(alignment)



