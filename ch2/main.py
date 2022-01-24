### Eli Blaney ###
###  BIO 501   ###
### Chapter 2  ###
### 2022-01-14 ###
print("Chapter 2 Exercises\n")

sequence = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"

# Output some basic data about the sequence for my own reference.
print("Sequence:", sequence)
print("Sequence length:", len(sequence))

### Exercise 2.1 ###
print("\nExercise 2.1\n")

# Determine the number of A and T bases in the sequence separately
num_a = sequence.count("A")
num_t = sequence.count("T")

# Add them together and find the percentage of AT bases in the sequence
num_at = num_a + num_t
at_percent = num_at / len(sequence) * 100

print("Total A:", num_a)
print("Total T:", num_t)
print("Total AT:", num_t)
# .format() is a method that embeds information into strings
# {:.2f} means to round the given float to 2 decimal places
print("AT%: {:.2f}%".format(at_percent))

### Exercise 2.2 ###
print("\nExercise 2.2\n")

# Find the complement of the sequence by replacing each base with its pair
# Because you want to avoid overwriting important data, turn them into lowercase so they
# don't conflict with one another.
# After that's done, turn the entire thing back to uppercase.
# Because each method returns a new string, I can just chain them together and make things really simple/easy.
complement = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").upper()
print("Complement:", complement)

### Exercise 2.3 ###
print("\nExercise 2.3\n")

sequence = "ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT"

# This is the sub-sequence that the restriction enzyme will bind to
match_seq = "GAATTC"
# This number represents which index in the match_seq that the restriction enzyme will actually cut at
cut_offset = 1

# Find where the restriction enzyme binds (only the first binding site)
seq_index = sequence.find(match_seq)

# Generate 2 fragments based on where it cuts by splitting the string
fragment_1 = sequence[:seq_index + cut_offset]
fragment_2 = sequence[seq_index + cut_offset:]

print("Fragment 1:", fragment_1)
print("Fragment 1 length:", len(fragment_1))
print("Fragment 2:", fragment_2)
print("Fragment 2 length:", len(fragment_2))

### Exercise 2.4a ###
print("\nExercise 2.4\n")

sequence = "ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"

# These are the indices that represent where the first exon ends and the second exon starts
# The first exon is from the beginning to base 63 (index 62)
# The intron is in the middle
# The second exon is from base 91 (index 90) to the end
exon_1_end = 62
exon_2_start = 90

# Generate the two exons as strings absed on their indices
exon_1 = sequence[:exon_1_end]
exon_2 = sequence[exon_2_start:]

print("Exon 1:", exon_1)
print("Exon 1 length:", len(exon_1))
print("Exon 2:", exon_2)
print("Exon 2 length:", len(exon_2))

### Exercise 2.4b ###
print("\nExercise 2.4b\n")

# Determine the total length of the exons combined (i.e. coding region)
coding_len = len(exon_1) + len(exon_2)

# Find the percent of the sequence that is coding
total_len = len(sequence)
coding_percent = coding_len / total_len * 100

print("Coding %: {:.2f}%".format(coding_percent))

### Exercise 2.4c ###
print("\nExercise 2.4c\n")

# Find the intron by splitting the sequence at the exon boundaries
intron = sequence[exon_1_end:exon_2_start]

# Build the sequence again, but this time making the intron lowercase to differentiate it.
new_seq = exon_1 + intron.lower() + exon_2

print("New sequence:", new_seq)