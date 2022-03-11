# The day has come -- I'm importing numpy for pulling off the On Your
# Own Project extra exercise #2. Only used it for that one in particular,
# because matrix multiplication sucks to do by hand.
import numpy as np
import math

# First, we need to find out what input files contain the 
# sequences that we are interested in. We'll accept
# them as paths from the user.
inputs = []
user_input = "ignored"
# Until the user pressees enter, we can continue to accept
# files, because "\n" is falsey, but file paths are truthy
while user_input:
    user_input = input("Enter input file paths (or enter to stop): ")
    if user_input:
        inputs.append(user_input.strip())

# Read all the sequences in paired order from the files
seqs = []
for filename in inputs:
    file = open(filename)
    # Ignoring FASTA comment lines, read every line in the file and parse it to uppercase
    seqs.extend([seq.rstrip().upper() for seq in file if seq and seq[0] != '>'])
    file.close()

# Grab the used amino acid symbols by parsing alphabetic characters from every sequence
symbols = [char for char in set(''.join(seqs)) if char.isalpha()]
# Convert the sequences to tuple pairs by zipping them with themselves at an offset of one
seqs = list(zip(seqs[::2], seqs[1::2]))

def get_score(paired_seqs, a1, a2):
    """Calculate the log-odds score for two amino acids.

    Keyword arguments:
    paired_seqs -- The list of tuple-paired sequences to train with
    a1 -- The first amino acid
    a2 -- The second amino acid
    """
    # Parameters to gather about the paired sequences and their
    # relationship with the amino acids in question.
    num_i = 0
    num_j = 0
    num_ij_aligned = 1 # Pseudocount
    num_ungapped_aligned_positions = 0
    for seq1, seq2 in paired_seqs:
        # Length is same for seq1 and seq2 since they are aligned
        length = len(seq1)
        # Count frequency of each amino acid in both sequences
        num_i += seq1.count(a1) + seq2.count(a1)
        num_j += seq1.count(a2) + seq2.count(a2)

        # Count ungapped aligned amino acid positions and the
        # frequency of both amino acids being aligned together.
        for k in range(length):
            if seq1[k] == seq2[k]:
                num_ij_aligned += 1
            if seq1[k] != '-' and seq2[k] != '-':
                num_ungapped_aligned_positions += 1

    # Number of total ungapped positions is considers both sequences, so double
    num_ungapped_positions = num_ungapped_aligned_positions * 2

    # Parameters to calculate the log-odds score
    # This formula is described on pages 98-99 of the Exploring Bioinformatics book
    qij = num_ij_aligned / num_ungapped_aligned_positions
    pi = num_i / num_ungapped_positions
    pj = num_j / num_ungapped_positions
    eij = 2 * pi * pj
    if a1 == a2:
        eij = pi ** 2
    sij = math.log2(qij / eij)

    return sij

# Convert the list of symbols into a matrix by filling in the scores for
# each paired combination of amino acid symbol.
matrix = [[get_score(seqs, a1, a2) for a2 in symbols] for a1 in symbols]

coli_power = input("\nCOLI value [1]: ")
if coli_power:
    # Convert to a numpy array matrix for easy matrix multiplication and then back to a list
    matrix = np.linalg.matrix_power(np.array(matrix), int(coli_power.strip())).tolist()
        
# Handy matrix printing function taken from my other scripts
print_matrix=lambda m:print('-'*(7*len(m[0])+3),*['| '+' '.join(map(lambda x:' '*(6-len(x))+x,map(lambda x:'{:.1f}'.format(x),r)))+' |'for r in m],'-'*(7*len(m[0])+3),sep='\n')if m else print('(empty matrix)')

# Print the generated information for the user to review
print('\nSymbols: ' + ', '.join(symbols))
print_matrix(matrix)

# Output matrix to a file that resembles a FASTA-like format
name = input('\nEnter name of matrix to save: ')
if name:
    filename = name.rstrip().lower().replace(' ', '_') + '.txt'
    file = open(filename, 'w')
    file.write('>' + name + '\n')
    file.write(','.join(symbols) + '\n')
    [file.write(','.join(map(str, row)) + '\n') for row in matrix]
    file.close()
    print('Saved to ' + filename)