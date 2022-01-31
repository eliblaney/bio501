import re

# Prepare input file for reading and output file for writing
data = open('Student_DATA.txt')
output = open('Student_DATA_cleaned.txt', 'w')

# Here I'm using a dict to speed up coding and optimize finding sequences
# using regular expressions. We're gonna use it for all the sequences,
# so I'm writing it just once outside the loop.
colors = { 'ATTAAGC': 'blue', 'ATTAGGC': 'brown', 'ATGGAGC': 'black', 'TTTAAGC': 'white' }
# Build a regular expression that allows us to match ANY of the keys
# in our dict from any given sequence.
pattern = '(' + '|'.join(colors.keys()) + ')'

# Go through all of the students' sequences
for line in data:
    # Split up the data to 2 variables
    name, seq = line.rstrip().split(",")
    # Make sequence uppercase, remove primers, and remove invalid characters
    seq = re.sub(r'[^ATCG]', '', seq[5:-5].upper())
    # Output the cleaned sequence to the requested file
    output.write(name + ',' + seq + '\n')
    # Find any of the color sequences from our dict
    m = re.search(pattern, seq)
    if m:
        # Figure out the colony color and print it
        print(name, 'has', colors.get(m.group(1), 'unknown'), 'colored colonies')
    else:
        # Missing a sequence for color? Say that.
        print(name, 'has no sequences corresponding to colony color.')
    # Finally, find and print GC content
    gc_content = (seq.count('G') + seq.count('C')) * 100.0 / len(seq)
    print('{} sequence has a GC content of {:.2f}%\n'.format(name, gc_content))

# We don't need these files anymore.
data.close()
output.close()
