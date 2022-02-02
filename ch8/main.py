### Eli Blaney ###
###  BIO 501   ###
### Chapter 8  ###
### 2022-02-01 ###
print('Chapter 8 Exercises')

import re

# As a challenge to myself and to one-up Gavril so I can poke fun at him,
# I decided to do each of the assigned exercises in one line. Obviously,
# the code is hard to read and a bit inefficient as a result.

### Bonus Exercise 8.0 ###
print('Bonus Exercise 8.0')

# This operates pretty much on a list of 3 tuples that outline the
# instructions for how to build the 3 desired dictionaries: pairing
# nucleotides with 3-letter amino acid codes, nucleotides with 1-
# letter amino acid codes, and 1-letter amino acid codes with 3-
# letter amino acid codes. The lambda functions dictate how to
# arrange the pairings for each of the three dictionaries. Next,
# each iteration opens and reads the raw gene codes and parses the
# text using regex to separate out each of the three components of
# the nucleotide, 1 letter code, and 3 letter codes. Then, the
# lambda functions assign keys and values for each of the three
# dictionaries, and those dictionaries are stored in another dictionary
# which simply names them so they are easy to access, and that dict
# is simultaneously assigned to the variable `ds` as well as printed.
print(ds:={n:{k(e):l(e) for e in re.findall('(.{3}).(.).(.{3}).?',open('genecode_raw1.txt').read())} for n,k,l in [('C1',lambda e:e[0],lambda e:e[1]),('C3',lambda e:e[0],lambda e:e[2]),('13',lambda e:e[1],lambda e:e[2])]})

### Exercise 8.1 ###
print('Exercise 8.1')

# Much simpler code, reusing the `ds` variable generated by the previous
# exercise. We accept a DNA sequence as input from the user, then make
# it uppercase and remove anything that's not a base pair. Next, we split
# up the DNA sequence to 3 letter substrings which represent codons. The
# codons are passed into the trinucleotide-1-letter-amino-acid dictionary
# created in the previous exercise and combined using the .join() method.
print('\nTranslated protein sequence:', ''.join([ds['C1'][codon] for codon in re.findall(r'.{3}', re.sub(r'[^ATCG]', '', input("\nEnter a DNA sequence: ").upper()))]))