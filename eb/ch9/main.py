# I've decided to use regex to parse the ORFs, rather than iteration
import re
import os
import sys
import time

def rc(seq):
    """Find the reverse complement of a sequence.

    Keyword arguments:
    seq -- A DNA sequence
    """
    return seq[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

def print_write(file, s):
    """Helper function to both print to stdout as well as write to a file.

    Keyword arguments:
    file -- The file to write to
    s -- The string to print and write
    """
    print(s.strip())
    file.write(s)

def find_orfs(genome, threshold=100, start_codon=r"ATG", stop_codon=r"(?:TAG|TGA|TAA)"):
    """Finds and sorts all the ORFs in a genome on the positive and negative strands.

    Keyword arguments:
    genome -- The genome to search
    threshold -- The minimum number of amino acids allowed in an ORF
    start_codon -- The pattern for a start codon. Defaults to r"ATG"
    stop_codon -- The pattern for a stop codon. Defaults to r"(?:TAG|TGA|TAA)"
    """
    genome_length = len(genome)

    # Build the regex pattern for an ORF
    # It searches for "ATG", then lazily searches for at least the
    # threshold number of codons (minus 1 for the start codon), then
    # any of the provided stop codons. The lazy search ensures that the
    # ORF terminates at the first stop codon that is found in frame.
    orf_pattern = start_codon + "(?:(?!" + stop_codon + ").{3}){" + str(threshold - 1) + ",}?" + stop_codon

    orfs = []
    avg = 0

    # Check the positive strand
    for m in re.finditer(orf_pattern, genome):
        start = m.start()
        end = m.end()
        length = end - start
        avg += length

        orf = genome[start:end]
        orfs.append((start + 1, '+', length, orf))

    # Check the negative (reverse complement) strand
    rc_genome = rc(genome)
    for m in re.finditer(orf_pattern, rc_genome):
        start = m.start()
        end = m.end()
        length = end - start
        avg += length

        orf = rc_genome[start:end]
        # Position corresponds with the positive strand, starting
        # with the final nucleotide in the reverse complement ORF.
        orfs.append((genome_length - end, '-', length, orf))

    # Calculate total, average, and sort the ORFs by nucleotide position
    total_orfs = len(orfs)
    if total_orfs > 0:
        avg /= total_orfs
        orfs = sorted(orfs, key=lambda o: o[0])

    return orfs, total_orfs, avg

print("To use [default values], just press enter rather than typing a value.\n")

# Grab the input file name for the genome from the user
filename = input('Genome file [W3110.txt]: ').strip() or "W3110.txt"

if not os.path.exists(filename):
    print("\nError: The provided file \"{}\" does not exist.".format(filename))
    sys.exit(1)

# Grab the amino acid minimum threshold from the user
threshold = 100

try:
    threshold_str = input('Minimum AAs [100]: ')
    if threshold_str:
        threshold = int(threshold_str)
        if threshold < 1:
            raise Exception()
except:
    print("\nError: The threshold value must be a positive integer.")
    sys.exit(1)

fasta_filename = 'output.fasta'
orfs_filename = 'output.txt'

# Read the input file and extract the genome
file = open(filename)
lines = list(file)
genome = ''.join(map(lambda x: x.strip(), lines[1:]))
file.close()

# Extra blank line to separate parameters from final output
print('')

# Get all the ORFs in sorted order and measure how long it takes
start_time = time.time()
orfs, total_orfs, avg = find_orfs(genome, threshold=threshold)
end_time = time.time()

print("Finished in {:.2f} seconds.".format(end_time - start_time))

if not total_orfs:
    print('The provided genomes does not contain any ORFs with {} minimum amino acids'.format(threshold))
else:
    # Now we need to output all the data
    # This one is super easy, just print out the ORF object in FASTA format
    with open(fasta_filename, 'w') as file:
        for o in orfs:
            file.write('> {} {} {}\n{}\n'.format(*o))

    # To output the overview details, we'll pad each column to 10 characters
    # We also want to print out the first 5 entires to the screen, so
    # we take advantage of the earlier defined print_write() function.
    pad = lambda l: ''.join([str(x) + ' ' * (10-len(str(x))) for x in l])
    with open(orfs_filename, 'w') as file:
        print_write(file, 'The provided genome contains {} ORFs with an average length of {:.0f} nucleotides.\n\n'.format(total_orfs, avg))
        print_write(file, pad(['Location','Strand','Length']))
        for o in orfs[:5]:
            print_write(file, '\n' + pad(o[:3]))
        # For the remaining items, only write them to the file.
        for o in orfs[5:]:
            file.write('\n' + pad(o[:3]))

    # Add ellipses and explain which files were generated
    print(pad(['...']*3) + '\n\nFull output written to', orfs_filename, 'and', fasta_filename)