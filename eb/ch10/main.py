# Several imports for fun shenanigans
import math
import re
import random
import sys

try:
    from matplotlib import pyplot as plt
except ImportError as e:
    print("[WARNING] You do not have matplotlib installed, so")
    print("[WARNING] this program will not display any graphs.")
    print("[WARNING]")
    print("[WARNING] If you have pip, you can install it with the following command:")
    print("[WARNING]        pip install matplotlib")

try:
    from statsmodels.stats import weightstats
except ImportError as e:
    print("[ERROR] You must install the statsmodels package to run this program.")
    print("[ERROR]")
    print("[ERROR] If you have pip, you can install it with the following command:")
    print("[ERROR]        pip install statsmodels")
    sys.exit(1)

# Standard values
amino_acids = "ARNDCQEGHILKMFPSTWYV*"
codon_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
            
filename = input("Filename [output.fasta]: ") or "output.fasta"

# Read output from Chapter 9 data
file = open(filename)
orfs = []
comment = None
for line in file:
    if line.startswith(">"):
        line = line.strip()[2:].split(" ")
        comment = [int(line[0]), line[1], int(line[2])]
    else:
        orfs.append((*comment, line.strip()))
file.close()

print("\n--------------------------------\nExercise 10.1\n")

def count_overlaps(orfs):
    """Count the number of overlapping regions in a list of ORFs.

    Keyword arguments:
    orfs -- The list of ORF tuples, including start pos, strand, length, and ORF contents
    """
    overlaps = 0
    last_pos_end = -1
    last_neg_end = -1
    for start, strand, length, orf in orfs:
        if strand == "+":
            last_pos_end = start + length
            if start < last_neg_end:
                overlaps += 1
        else:
            last_neg_end = start + length
            if start < last_pos_end:
                overlaps += 1
    return overlaps

num_overlaps = count_overlaps(orfs)

print("Number of overlapped regions: {}".format(num_overlaps))

print("\n--------------------------------\nExercise 10.2\n")

def translate(seq):
    """Translate a DNA sequence to amino acids.

    Keyword arguments:
    seq -- The DNA sequence to translate
    """
    return ''.join([codon_table[codon] for codon in re.findall('.{3}', re.sub('[^ATCG]', '', seq))])

def freq(seqs, aa):
    """Determine the frequency of a single amino acid in a list of sequences.

    Keyword arguments:
    seqs -- A list of amino acid sequences
    aa -- The single-letter code for the amino acid to search
    """
    length = 0
    count = 0
    for seq in seqs:
        length += len(seq)
        count += seq.count(aa)
    return count, count / length

# Get a list of amino acid sequences and find the frequencies for all amino acids
orf_seqs = [o[3] for o in orfs]
seqs = [translate(o) for o in orf_seqs]
freqs = {aa: freq(seqs, aa) for aa in amino_acids}

print("Freqency of each amino acid:")
print('AA\tCount\t%\tDiff from 5%')
[print('{}\t{}\t{:.2f}%\t{:.2f}'.format(aa, f[0], f[1] * 100, f[1] * 100 - 5)) for aa, f in freqs.items()]
total = sum(map(lambda f: f[1], freqs.values())) * 100
print("========================================")
print('Total: {:.0f}%\tMean diff from 5%: {:.2f}'.format(total, total / 20 - 5))

print("\n--------------------------------\nExercise 10.3\n")

def generate_orf(length, bases='ATCG', weights=[1, 1, 1, 1]):
    """Generate a random ORF nucleotide sequence.

    Keyword arguments:
    length -- The length of the ORF to generate
    bases -- The bases available to use. Defaults to 'ATCG'.
    weights -- The probability weights for each base. Defaults to [1, 1, 1, 1].
    """
    return ''.join(random.choices(bases, weights=weights, k=length))

def gc(orfs):
    """Convert ORFs to a list of GC percentages.

    Keyword arguments:
    orfs -- The list of ORF sequences to convert
    """
    gcs = []
    for o in orfs:
        g = o.count('G')
        c = o.count('C')
        gcs.append((g + c) / len(o))
    return gcs

# Create a null distribution from the randomly generated ORFs
null_orfs = [generate_orf(len(o)) for o in orf_seqs]
# Then convert both ORF lists to GC% distributions
null_gc = gc(null_orfs)
orf_gc = gc(orf_seqs)

print('Null ORF GC%: {:.4f}'.format(sum(null_gc) / len(null_gc) * 100))
print('Genic ORF GC%: {:.4f}'.format(sum(orf_gc) / len(orf_gc) * 100))

# Run a one-tailed t-test
t, p, _ = weightstats.ttest_ind(orf_gc, null_gc)
p /= 2
print("t-statistic: {:.4f}".format(t))
print("One-tailed p-value: {:.4f}".format(p))
print("             ln(p): {:.4f}\n".format(math.log(p)))

if p < 0.05:
    print("The p-value is less than 0.05, so the null hypothesis is rejected. Genic GC% is significantly higher than inter-genic areas of the genome.")
else:
    print("The p-value is at least 0.05, so the null hypothesis cannot be rejected. Genic GC% is NOT significantly higher than inter-genic areas of the genome.")

# Optional graphing if matplotlib is installed and the user wants it
if 'matplotlib' in sys.modules and 'n' != input('\nShow intra/inter-genic distributions [y]? '):
    plt.hist(null_gc, 250, fc=(0, 0, 1, 0.8))
    plt.hist(orf_gc, 250, fc=(1, 0.5, 0, 0.8))
    plt.show()

print("\n--------------------------------\nExercise 10.4\n")

max_dist = int(input('Operon max distance [25]: ') or 25)
threshold = int(input('Operon min streak [2]: ') or 2)

# Count streaks of nearby ORFs on both strands to determine what operons exist
streaks = []
for strand in ['+', '-']:
    last_end = 0
    streak = 1
    for o in [o for o in orfs if o[1] == strand]:
        if o[0] - last_end < max_dist:
            streak += 1
        else:
            streaks.append(streak)
            streak = 1
        last_end = o[0] + o[2]

# Determine operons based on if the ORF streaks are above the user-specified threshold
operons = [s for s in streaks if s >= threshold]
num_in_operons = sum(operons)

print('\nAverage streak: {:.4f} ORFs'.format(sum(streaks) / len(streaks)))
print('Number of predicted operons: {}'.format(len(operons)))
print('Num in operons: {}'.format(num_in_operons))
print('Fraction in operons: {:.2f}%\n'.format(num_in_operons / len(orfs) * 100))

# Run a one-tailed z-test using the streaks distribution and
# a null distribution filled with ones
z, p = weightstats.ztest(streaks, value=1)
p /= 2

print('Z-statistic: {:.2f}'.format(t))
print('p-value: {:.2f}'.format(p))
print("  ln(p): {:.2f}\n".format(math.log(p)))

if p < 0.05:
    print("The p-value is less than 0.05, so the null hypothesis is rejected. There is significant evidence to say that many E. coli genes are in operons.")
else:
    print("The p-value is at least 0.05, so the null hypothesis cannot be rejected. There is NOT significant evidence to say that many E. coli genes are in operons.")
