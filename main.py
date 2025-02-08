from Bio import SeqIO
from Bio.Align import PairwiseAligner

# load the reference FASTA sequence file
reference_record = SeqIO.read("/Users/albertseo/Development/covid19-mutations/sequences/EPI_ISL_402124-ORF1ab.fasta", "fasta")
reference_seq = reference_record.seq

# load the mutant sequences (can accept more than one variant sequence)
mutant_records = list(SeqIO.parse("/Users/albertseo/Development/covid19-mutations/sequences/MZ571142.1.fasta", "fasta"))

# align the sequences to detect differences using PairwiseAligner
# PairwiseAligner provides more flexibility and better control over alignment parameters.
aligner = PairwiseAligner() # establish the object
aligner.mode = "global"  # global (Needleman-Wunsch) alignment
aligner.match_score = 1 # a score that indicates how well a read nucleotide matches a reference nucleotide
aligner.mismatch_score = -1 # the negative penalty assigned when two characters from different sequences are aligned that are not identical
aligner.open_gap_score = -2 # the penalty value assigned to the act of initiating a gap (inserting a space) when aligning two sequences
aligner.extend_gap_score = -0.5 # the penalty assigned to extending a gap by one base in a sequence alignment

for mutant_record in mutant_records:
    mutant_seq = mutant_record.seq
    alignments = aligner.align(reference_seq, mutant_seq)

    print(f"\nAlignment for {mutant_record.id}:")
    print(alignments[0])  # Display best alignment

# detect mutations by detecting single nucleotide polymorphisms (SNPs), insertions, and deletions.
# Compare each aligned query sequence with the reference
for mutant_record in mutant_records:
    mutant_seq = mutant_record.seq
    mutations = [] # extract and initialize a mutant list

    # compare each nucleotide position while iterating over each sequence
    for i, (ref_nuc, mutant_nuc) in enumerate(zip(reference_seq, mutant_seq)):
        if ref_nuc != mutant_nuc:  # if the nucleotides are different then it is a mutant
            # classify the mutant if it is detected
            if ref_nuc == "-":
                mutation_type = "Insertion"
            elif mutant_nuc == "-":
                mutation_type = "Deletion"
            else:
                mutation_type = "SNP"
            mutations.append(f"Pos {i+1}: {ref_nuc} → {mutant_nuc} ({mutation_type})") # store the mutant information

    print(f"\nMutations in {mutant_record.id}:")
    print("\n".join(mutations) if mutations else "No mutations found")