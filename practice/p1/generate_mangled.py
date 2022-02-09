import random

grape_promoters = open('grape_promoters.txt')
mangled_promoters = open('mangled_promoters.txt', 'w')

random_chars = "BDEFHIJKLMNOPQRSUVWXYZ"

for line in grape_promoters:
    gene, seq = line.rstrip().split(",")
    seq = ''.join(random.choice((str.upper, str.lower))(base) for base in seq)
    len_seq = len(seq)
    for _ in range(random.randrange(0, 100)):
        i = random.randrange(0, len_seq)
        ch = random.choice(random_chars)
        seq = seq[:i] + ch + seq[i:]
    mangled_promoters.write(gene + "," + seq + "\n")

grape_promoters.close()
mangled_promoters.close()
