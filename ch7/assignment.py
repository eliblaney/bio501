import re

grape_promoters = open('mangled_promoters.txt')

def get_num_at(dna):
    # dna = re.sub(r"[^ATCG]", "", dna.upper())
    # a = dna.count('A')
    # t = dna.count('T')
    # return a + t
    return len(re.sub(r"[^AT]", "", dna.upper()))

gene_most_at = ""
most_at = -1
for line in grape_promoters:
    gene, seq = line.rstrip().split(",")
    at = get_num_at(seq)
    if at > most_at:
        gene_most_at = gene
        most_at = at

print(gene_most_at, most_at)

grape_promoters.close()