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

accession_names = ['xkn59438', 'yhdck2', 'eihd39d9', 'chdsye847', 'hedle3455', 'xjhd53e', '45da', 'de37dp']
patterns = [r'5', r'[de]', r'de', r'd\we', r'(de|ed)', r'^[xy]', r'^[xy].*e$', r'\d{3,}', r'd[arp]$']

for pattern in patterns:
    print('\nPattern:', pattern)
    for name in accession_names:
        if re.search(pattern, name):
            print(name, end=' ')

### Exercise 7.2 ###
print('\n\nExercise 7.2\n')

dna_file = open('dna.txt')
dna = dna_file.read().rstrip().upper()
dna_file.close()

def cut(seq, enzyme, cut_site, complement=True):
    """Determine the length of fragments created when a DNA sequence is cut with a restriction enzyme.

    Keyword arguments:
    seq -- The DNA sequence to cut
    enzyme -- The enzyme's regular expression recognition site
    cut_site -- The index of the enzyme's cut site relative to the regex
    complement -- Whether to consider the complement DNA fragment's potential for being cut. Defaults to True.

    """
    length = len(seq)
    cuts = []
    matches = re.finditer(enzyme, seq)
    for m in matches:
        cut_index = m.start() + cut_site
        cuts.append(cut_index)
    if complement:
        complement = seq[::-1].replace('A', 't').replace('T', 'a').replace('G', 'g').replace('C', 'c').upper()
        matches = re.finditer(enzyme, complement)
        for m in matches:
            cut_index = length - (m.start() + cut_site) + 1
            cuts.append(cut_index)

    fragment_sizes = []
    last = 0
    cuts.sort()
    for cut in cuts:
        fragment_sizes.append(len(seq[last:cut]))
        last = cut
    fragment_sizes.append(len(seq[last:]))

    return fragment_sizes

AbcI = r'A.TAAT'
AbcI_cut = 3
AbcII = r'GC[AG][AT]TG'
AbcII_cut = 4

print('AbcI makes fragments with lengths', cut(dna, AbcI, AbcI_cut, False))
print('AbcI makes fragments with lengths', cut(dna, AbcI, AbcI_cut), 'when considering complement')
print('AbcI makes fragments with lengths', cut(dna, AbcII, AbcII_cut, False))
print('AbcI makes fragments with lengths', cut(dna, AbcII, AbcII_cut), 'when considering complement')
