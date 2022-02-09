### Eli Blaney ###
###  BIO 501   ###
### Chapter 9  ###
### 2022-02-03 ###
print('Chapter 9 Exercises')

# Import our operating system and system libraries
import os
import sys

### Exercise 9.1 ###
print('\nExercise 9.1\n')

# These parameters can be changed but are the default
# folders to store the input and output DNA sequence files
input_files_dir = 'input_files'
output_files_dir = 'output_files'
# This is the system file separator, default '/' for *nix
# and macOS systems, and '\' for Windows systems
file_separator = '/'
# The factor of which to divide our sequence to find out which
# bin the sequence should go into.
bin_factor = 100

# Iterate through all of the files in the input directory
for filename in os.listdir(input_files_dir):
    # Build the path to each file by combining the
    # input folder path with the file name
    filepath = input_files_dir + file_separator + filename
    file = open(filepath)

    # Sequence counter so we can keep track of the sequences easily.
    # Not strictly required for this assignment, but I like it.
    i = 0
    # Here, we iterate through all the sequences in each file to
    # sort each of them into the correct bin
    for seq in file:
        seq = seq.rstrip()
        # Increment our sequence counter
        i += 1
        # Find the length and use INTEGER DIVISION to find out
        # which bin the DNA sequence goes to. Integer division
        # involves normal division, and then flooring the result
        # to the great integer below the quotient.
        # This corresponds with the bin we want to insert the
        # sequence into.
        length = len(seq)
        bin = length // bin_factor
        # Build the path to the output folder
        output_folder = output_files_dir + file_separator + str(bin)
        # Create the bin folder if it isn't already there
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        # Build the output file path by adding on the original
        # file name to the output folder path. The sequence number
        # of the original file is also added on with an underscore.
        # That way, we can easily relocate certain sequences.
        output_path = output_folder + file_separator + str(i) + '_' + filename
        # Open, write, and close the output file.
        output_file = open(output_path, 'w')
        output_file.write(seq)
        output_file.close()

    print("[{}] Sorted {} sequences.".format(filename, i))
    file.close()

### Exercise 9.2 ###
print('\nExercise 9.2\n')

# Accept arguments from the command line, OR
# we can also use default arguments if they
# are not provided. There are libraries to do
# this more effectively, but this strategy
# works completely well for this assignment.
if len(sys.argv) >= 3:
    kmer_length = int(sys.argv[1])
    cutoff = int(sys.argv[2])
if len(sys.argv) < 3:
    cutoff = 900
    print("Using default cutoff of", cutoff)
if len(sys.argv) < 2:
    kmer_length = 3
    print("Using default kmer_length of", kmer_length, '\n')

# Frequency dictionary for each k-mer
freqs = {}
# Iterate through all the files in our input folder
for filename in os.listdir(input_files_dir):
    # Build the file path like before
    filepath = input_files_dir + file_separator + filename
    file = open(filepath)
    # For every sequence in the file, find every k-mer of
    # the given length
    for seq in file:
        seq = seq.rstrip()
        # This range trick allows us to easily iterate through
        # and grab substrings that correspond with each k-mer
        # without thinking too hard.
        for i in range(kmer_length, len(seq)):
            kmer = seq[i - kmer_length:i]
            # Insert it into the dictionary or increment
            # its frequency if it already exists.
            if kmer in freqs:
                freqs[kmer] = freqs[kmer] + 1
            else:
                freqs[kmer] = 1

    file.close()

# Now, we can filter the k-mer frequencies to
# only print out the k-mers that are above
# our cutoff variable.
print("{}-mers with frequency more than {}:".format(kmer_length, cutoff))
for kmer in freqs.keys():
    freq = freqs[kmer]
    if freq > cutoff:
        print("{}: {}".format(kmer, freq))

