# Dictionary of amino acids and their hydrophobicities
amino_acids = {
    'D': -3.5,
    'K': -3.9,
    'H': -3.2,
    'T': -0.7,
    'V':  4.2,
    'F':  2.8,
    'M':  1.9,
    'Y': -1.3,
    'N': -3.5,
    'L':  3.8,
    'E': -3.5,
    'R': -4.5,
    'W': -0.9,
    'Q': -3.5,
    'I':  4.5,
    'C':  2.5,
    'A':  1.8,
    'S': -0.8,
    'G': -0.4,
    'P': -1.6
}

thresholds = [3, 0, -1]
str_thresholds = input('Thresholds [3,0,-1]: ')
if str_thresholds:
    thresholds = list(map(float, str_thresholds.split(',')))

def get_bucket(v, thresholds):
    """Get the bucket number that a value would fall into.

    Buckets are "categories" that values may fall into, depending
    on which range of threshold values they fall into.

    Keyword arguments:
    v -- The value to find a bucket number for.
    thresholds -- A descending array of thresholds to compare ranges.
    """
    # i is current the bucket number
    i = 0
    for t in thresholds:
        # If the value is greater than any of the threshold values,
        # then the bucket has been found
        if v > t:
            return i
        i += 1
    # If lower than all thresholds, it is the last bucket
    return i

# Convert the amino acid hydrophobicities to a dictionary of buckets
buckets = { k: get_bucket(v, thresholds) for k, v in amino_acids.items() }

def get_score(a1, b1, a2, b2):
    """Get the substutition score of two amino acids.

    Keyword arguments:
    a1 -- The first amino acid letter
    b1 -- The first amino acid bucket
    a2 -- The second amino acid letter
    b2 -- The second amino acid bucket
    """
    # Match = 2
    if a1 == a2:
        return 2
    # Same bucket = 1
    elif b1 == b2:
        return 1
    # Similar hydrophobicity or hydrophilicity = 0.5
    elif (b1 > 1 and b2 > 1) or (b1 <= 1 and b2 <= 1):
        return 0.5
    # Mismatch = 0
    return 0

# Build the matrix
matrix = [[get_score(a1, b1, a2, b2) for a2, b2 in buckets.items()] for a1, b1 in buckets.items()]

print_matrix=lambda m:print('-'*(5*len(m[0])+3),*['| '+' '.join(map(lambda x:' '*(4-len(x))+x,map(lambda x:'{:.1f}'.format(x),r)))+' |'for r in m],'-'*(5*len(m[0])+3),sep='\n')if m else print('(empty matrix)')

print(list(buckets.keys()))
print_matrix(matrix)

name = input('\nEnter name of matrix to save: ')
if name:
    filename = name.rstrip().lower().replace(' ', '_') + '.txt'
    file = open(filename, 'w')
    file.write('>' + name + '\n')
    file.write(','.join(buckets.keys()) + '\n')
    [file.write(','.join(map(str, row)) + '\n') for row in matrix]
    file.close()
    print('Saved to ' + filename)