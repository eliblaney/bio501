### Eli Blaney ###
###  BIO 501   ###
### Chapter 3  ###
### 2022-01-19 ###
print("Chapter 3 Exercises\n")

### Exercise 3.1 ###
print("\nExercise 3.1\n")

# Open the genomic_DNA.txt file to get the sequence
genomic_dna = open("genomic_DNA.txt")
# Store the sequence from the file and remove trailing whitespace
sequence = genomic_dna.read().rstrip()
# Don't need the file anymore, close it
genomic_dna.close()

print("Sequence:", sequence)

# These are the indices that represent where the first exon ends and the second exon starts
# The first exon is from the beginning to base 63 (index 62)
# The intron is in the middle
# The second exon is from base 91 (index 90) to the end
exon_1_end = 62
exon_2_start = 90

# Generate the two exons as strings absed on their indices
exon_1 = sequence[:exon_1_end]
exon_2 = sequence[exon_2_start:]

# Build the sequence with only the exons included
exons = exon_1 + exon_2

print("Exons:", exons)

# Find the intron by splitting the sequence at the exon boundaries
intron = sequence[exon_1_end:exon_2_start]

print("Intron:", intron)

# Open the new file exons.txt for writing to print out the exons
exons_file = open("exons_3.1.txt", "w")
# Write out the combined exon sequence and close
exons_file.write(exons)
exons_file.close()

print("Wrote exons to exons_3.1.txt")

# Open the new file intron.txt for writing to print out the intron
intron_file = open("intron_3.1.txt", "w")
# Write and close like before
intron_file.write(intron)
intron_file.close()

print("Wrote intron to intron_3.1.txt")

### Exercise 3.2 ###
print("\nExercise 3.2\n")

# These are the sequences from the textbook stored as a list of tuples
# The first tuple item contains the header name, while the second
# contains the sequence itself
sequences = [
    ("ABC123", "ATCGTACGATCGATCGATCGCTAGACGTATCG"),
    ("DEF456", "actgatcgacgatcgatcgatcacgact"),
    ("HIJ789", "ACTGAC-ACTGT--ACTGTA----CATGTG")
]

# Open up our output file for the FASTA sequences
output = open("output_3.2.fasta", "w")

# Iterate over each of the tuple pairs to write all of the sequences
for k, v in sequences:
    # Write out the header of the sequence and go to a new line
    output.write(">" + k + "\n")
    # Write out the sequence itself after making it uppercase and removing invalid -'s
    output.write(v.upper().replace("-", "") + "\n")
    # Repeat the process for each sequence

# Once finished, close the output file
output.close()

print("Wrote sequences to output_3.2.fasta")

### Exercise 3.3 ###
print("\nExercise 3.3\n")

# Like before, iterate through the sequences
for k, v in sequences:
    # But this time, let each sequence open its own file
    output = open(k + "_3.3.fasta", "w")
    # Write out the header and sequence to its dedicated file
    output.write(">" + k + "\n")
    output.write(v.upper().replace("-", "") + "\n")
    # Close this sequence's file
    output.close()
    
    # Now we can print a success message and move on to the next sequence
    print("Wrote sequence to " + k + "_3.3.fasta")