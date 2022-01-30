import re

grape_promoters = open('mangled_promoters.txt')

target = "GSVIVT01020437001_15026"

def get_len_longest_at(dna):
    dna = re.sub(r"[^ATCG]", "", dna.upper())
    return len(re.findall(r"[AT]+", dna))

dna = ""
for line in grape_promoters:
    gene, seq = line.rstrip().split(",")
    if gene == target:
    at = get_len_longest_at(seq)
    if at > longest_at:
        gene_most_at = gene
        longest_at = at

print(gene_most_at, longest_at)

grape_promoters.close()