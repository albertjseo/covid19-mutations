from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd

# load the sequence files
gb_ref_file = "/Users/albertseo/Development/covid19-mutations/ncbi_dataset/data/GenBank/NC_045512.2.gb"
gb_var_file = "/Users/albertseo/Development/covid19-mutations/ncbi_dataset/data/GenBank/MZ571142.1.gb"

# parse the files but do not add .seq so it can be used for annotations as well
ref_seq = SeqIO.read(gb_ref_file, "genbank")
var_seq = SeqIO.read(gb_var_file, "genbank")

ref_seq1 = ref_seq.seq
var_seq1 = var_seq.seq

# checkpoint
print("Reference Seq:", ref_seq1)
print("Variant Seq:", var_seq1)


# perform alignment using aligner object
# PairwiseAligner provides more flexibility and better control over alignment parameters
aligner = PairwiseAligner() # establish the object
aligner.mode = "global"  # global (Needleman-Wunsch) alignment
aligner.match_score = 1 # a score that indicates how well a read nucleotide matches a reference nucleotide
aligner.mismatch_score = -1 # the negative penalty assigned when two characters from different sequences are aligned that are not identical
aligner.open_gap_score = -2 # the penalty value assigned to the act of initiating a gap (inserting a space) when aligning two sequences
aligner.extend_gap_score = -0.5 # the penalty assigned to extending a gap by one base in a sequence alignment

# access the best alignment
alignments = aligner.align(str(ref_seq), str(var_seq))[0]
print("Best Alignment:")
print(alignments)

# identify mutations
# pass aligned sequences to check for mutations
def identify_mutations(alignments):
    snps = []
    insertions = []
    deletions = []

    # Loop through the aligned sequences and identify SNPs, insertions, or deletions
    aligned_seq1 = alignments[0]
    aligned_seq2 = alignments[1]

    for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if base1 != base2:
            if base1 == '-' or base2 == '-': # Insertion in seq1
                insertions.append((i, base1, base2))
            elif base1 == '-' and base2 != '-':  # Deletion in seq2
                deletions.append((i, base1, base2))
            elif base1 != '-' and base2 != '-':
                snps.append((i, base1, base2))  # SNP (base difference)

    return snps, insertions, deletions

# Get SNPs, insertions, and deletions
snps, insertions, deletions = identify_mutations(alignments)

# parse gene annotations for later
def gene_annotation(record):
    genes = []

    for feature in record.features:
        if feature.type == "gene":
            gene_info = {
                "gene": feature.qualifiers.get("gene", ["Unknown"])[0],
                "start": feature.location.start,
                "end": feature.location.end,
            }
            genes.append(gene_info)
    return genes

# gene annotations
ref_seq_genes = gene_annotation(ref_seq)
var_seq_genes = gene_annotation(var_seq)

print("Reference Seq:", ref_seq_genes)
print("Variant Seq:", var_seq_genes)

# function that annotates mutations
def annotate_mutations(mutations, annotations):
    annotated_mutations = []
    for mutation in mutations:
        position, base1, base2 = mutation
        gene_found = False

        for gene in annotations:
            if gene["start"] <= position < gene["end"]:
                annotated_mutations.append((position, base1, base2, gene["gene"]))
                gene_found = True
                break
            if not gene_found:
                annotated_mutations.append((position, base1, base2, "Intergenic"))
    return annotated_mutations

annotated_snps = annotate_mutations(snps, ref_seq_genes)
annotated_insertions = annotate_mutations(insertions, ref_seq_genes)
annotated_deletions = annotate_mutations(deletions, ref_seq_genes)

# store the annotated mutations as a dataframe to use for visualizations later
snp_df = pd.DataFrame(annotated_snps, columns=["Position", "Seq1", "Seq2", "Gene"])
insertion_df = pd.DataFrame(annotated_insertions, columns=["Position", "Seq1", "Seq2", "Gene"])
deletion_df = pd.DataFrame(annotated_deletions, columns=["Position", "Seq1", "Seq2", "Gene"])
mutations_df = pd.concat([snp_df, insertion_df, deletion_df])

print("SNPs:")
print(snp_df)
print("\nInsertions:")
print(insertion_df)
print("\nDeletions:")
print(deletion_df)

