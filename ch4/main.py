### Eli Blaney ###
###  BIO 501   ###
### Chapter 4  ###
### 2022-01-23 ###
print("Chapter 4 Exercises\n")

### Exercise 4.1 ###
print("\nExercise 4.1\n")

# Open the input.txt file to get the sequences
input_dna = open("input.txt")
# Define the output_4.1.txt file object for trimmed DNA sequences
output_dna = open("output_4.1.txt", "w")

# The number of base pairs to strip from the beginning of each sequence
num_trim = 14

# Counter variable for nice printing
i = 1

# Iterate over every line (DNA sequence) in the input.txt
for sequence in input_dna:
    # Remove the first num_trim characters from the sequence
    # Here I use rstrip() as Dr. Cho suggested to remove the newlines and
    # give a more accurate sequence length
    trimmed = sequence[num_trim:].rstrip()

    # Print the original and trimmed sequence lengths
    # As an alternative to rstrip(), one could also simply subtract 1 from
    # the length, which could (potentially) be marginally faster,
    # but that optimization is not necessary for this kind of operation
    print("Original sequence", i, "length:", len(sequence))
    print(" Trimmed sequence", i, "length:", len(trimmed))

    # Write the trimmed sequence to the output_4.1.txt file
    # We also need to add the newline back since it was removed earlier
    output_dna.write(trimmed + "\n")
    # Increment the counter
    i += 1

# Close our files
input_dna.close()
output_dna.close()

### Exercise 4.2 ###
print("\nExercise 4.2\n")

# Open the genomic_dna.txt file containing the DNA sequences of interest
genomic_dna = open("genomic_dna.txt")
# Open the exons.txt file containing the indices that define the exons
exons = open("exons.txt")

# Parse and strip whitespace from the genomic DNA sequence
sequence = genomic_dna.read().rstrip()

# Parse the lines into a list of comma-separated indices
exon_indices = exons.read().rstrip().split("\n")

# We don't need the files anymore
exons.close()
genomic_dna.close()

# Open the output file for the parsed exons
output = open("output_4.2.txt", "w")

# Counter variable for nice printing
i = 1
for indices in exon_indices:
    # Split the comma-separated indices into distinct variables
    start_str, stop_str = indices.split(",")

    # The indices are strings, so we have to parse them as integers using int()
    # As Dr. Cho indicated, the provided start/stop positions are one-indexed and both inclusive
    # The start position needs to be decremented to be an inclusive index,
    # while the stop position needs no action to transform to an inclusive index
    start = int(start_str) - 1
    stop = int(stop_str)

    # We can use the integer indices to get the exon as a sequence substring
    exon = sequence[start:stop]

    # Print out the exon WITHOUT a newline.
    # This has the added bonus of automatically concatenating them
    output.write(exon)

    # Print out the parsed exon and increment the counter
    print("Exon", i, "=", exon)
    i += 1

# Write a final newline to the file
output.write("\n")

# Flush the written exon output to output_4.2.txt and close the file
output.close()

### Bonus Assignment ###
###   Exercise 4.3   ###
print("\nExercise 4.3 (Bonus Assignment)\n")

# Calculate the factorial of a number

# The number to calculate the factorial for
x = 5

# Define the output product
product = 1
# Loop over each number under x to multiply
for xi in range(x):
    # Multiply every number in the range 1:x, inclusive
    product = product * (xi + 1)

print("Factorial of", x, "is", product)