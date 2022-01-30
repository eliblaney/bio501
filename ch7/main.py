### Eli Blaney ###
###  BIO 501   ###
### Chapter 7  ###
### 2022-01-28 ###
print('Chapter 7 Exercises')

import re

### Exercise 7.1 ###
print('Exercise 7.1')

# Open the data file for reading
data = open('dna.txt')

# We can repeat the same code for all these names
accession_names = ['xkn59438', 'yhdck2', 'eihd39d9', 'chdsye847', 'hedle3455', 'xjhd53e', '45da', 'de37dp']
# We can also reuse the same code for all the regular expression patterns
patterns = [r'5', r'[de]', r'de', r'd\we', r'(de|ed)', r'^[xy]', r'^[xy].*e$', r'\d{3,}', r'd[arp]$']

# Iterate through all the patterns and accession names
for pattern in patterns:
    print('\nPattern:', pattern)
    for name in accession_names:
        # If we find the pattern in the name, print it out
        if re.search(pattern, name):
            # end is an optional keyword argument that lets us replace the newline character with
            # something else; in this case, a space, causing all the names to be on the same line.
            print(name, end=' ')

### Exercise 7.2 ###
print('\n\nExercise 7.2\n')

# Open the dna.txt file and extract the DNA sequence
dna_file = open('dna.txt')
# Make sure we strip away the newline and double check it's uppercase, so we can work
# with any given DNA sequence.
dna = dna_file.read().rstrip().upper()
# No need for the dna.txt file anymore
dna_file.close()

def cut(seqs, enzyme, cut_site, complement=True):
    """Determine the fragments created when a DNA sequence is cut with a restriction enzyme.

    Keyword arguments:
    seqs -- An array of DNA sequences to cut
    enzyme -- The enzyme's regular expression recognition site
    cut_site -- The index of the enzyme's cut site relative to the regex
    complement -- Whether to consider the complement DNA fragment's potential for being cut. Defaults to True.

    """
    # This will keep track of all the generated fragments, which we can return at the end
    total_fragments = []
    # Iterating through an array of sequences allows cuts to be chained consecutively,
    # passing the return value of this function as the seqs argument to its next call.
    for seq in seqs:
        # This keeps track of the indices where fragments will be split in our sequence
        cuts = []
        # Create an array of match objects to get our fragments
        matches = re.finditer(enzyme, seq)
        # Go through each match to find where it's cutting
        for m in matches:
            # Take the match's start index and add it to the enzyme's cut site
            # This gives us the index of our larger DNA sequence where the cut
            # will be made.
            cut_index = m.start() + cut_site
            # Store it for later.
            cuts.append(cut_index)
        # If we are considering complement, we can do a similar process, but reversed.
        if complement:
            # We need this to calculate the cut indices later.
            length = len(seq)
            # Create the complement by reversing the string and swapping out base pairs.
            # seq[::-1] is a special trick to reverse a string. It uses splicing to start
            # from the end of the string and work backwards to the beginning.
            complement = seq[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
            # Like before, find and iterate through the match objects
            matches = re.finditer(enzyme, complement)
            for m in matches:
                # Account for the fact that since the complement sequence is reversed,
                # we need to subtract the complement cut site index from the total
                # length and add one to find where the enzyme would cut in our normal
                # DNA sequence.
                cut_index = length - (m.start() + cut_site) + 1
                cuts.append(cut_index)

        # This array will keep track of the fragments for seq generated in this iteration.
        fragments = []
        # Tracker variable for cut indices as we iterate through them.
        last = 0
        # We need to sort all the indices since we are going to move forward through
        # the sequence.
        cuts.sort()
        # Now we can iterate through all of the cut indices, including complement if desired,
        # and generate the actual DNA fragments.
        for cut in cuts:
            # Slice from the last cut index to the current cut index.
            fragments.append(seq[last:cut])
            # Update the last cut index to this one.
            last = cut
        # We need to add on the final fragment that extends to the end of the sequence.
        fragments.append(seq[last:])
        # Finally, concat the fragments generated from this sequence to the larger
        # array containing all of the generated fragments from all sequences.
        total_fragments += fragments

    # After iterating through all of the sequences, we are left with all fragments generated
    # from the enzyme in all sequences.
    return total_fragments

def to_lengths(fragments):
    """Convert a list of sized elements to a list of the sizes of those elements."""
    fragment_lengths = []
    for f in fragments:
        # Straightforward: iterate through all fragments and add each
        # one's length to a newly generated list
        fragment_lengths.append(len(f))
    return fragment_lengths

# Now, the good part.
# Here are the two enzymes, represented as regular expressions.
# We also keep track of the cut index, which is an integer corresponding
# to the index in our regular expression that will be used for cutting
# a DNA sequence into two fragments. The index represents where the
# second fragment will begin.
AbcI = r'A.TAAT'
AbcI_cut = 3
AbcII = r'GC[AG][AT]TG'
AbcII_cut = 4

# To do a double digest, we simply call our cut function twice, chaining the output
# of the first call (cutting with AbcI) as our input to the second call (AbcII).
AbcI_fragments = cut([dna], AbcI, AbcI_cut, False)
AbcII_fragments = cut(AbcI_fragments, AbcII, AbcII_cut, False)
# We are given an actual array of fragments, which we can use for other purposes.
# However, here we will just transform the fragments into a list of the fragment
# sizes and print them out.
print('Double digest creates fragments with lengths', to_lengths(AbcII_fragments))

# Now we can do the same thing, but considering complement using the same function.
# The True argument for considering complement is optional as the function defaults to True,
# but I have included it here to be explicit about the difference between our function calls.
AbcI_fragments = cut([dna], AbcI, AbcI_cut, True)
AbcII_fragments = cut(AbcI_fragments, AbcII, AbcII_cut, True)
print('Double digest creates fragments with lengths', to_lengths(AbcII_fragments), "when considering complement")
