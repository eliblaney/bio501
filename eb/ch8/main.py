import random

def coverage_met(coverage, fold):
    i = 0
    met = True
    clen = len(coverage)
    while i < clen and met:
        met &= coverage[i] >= fold
        i += 1
    return met

def complement(seq):
    return seq[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

def build_fragments(seq, fmin, fmax, fold, error_rate=0.1):
    seqs = [seq, complement(seq)]
    length = len(seq)
    coverage = [0] * length * 2
    frags = []
    reads = 0
    while not frags or not coverage_met(coverage, fold):
        dna = random.choice(seqs)
        rlen = random.randrange(fmin, fmax + 1)
        rs = random.randrange(0, length + 1 - rlen)
        fragment = dna[rs:rs+rlen]
        if random.random() / (0.000001 + ((abs(rs - length/2)/length) ** 2)) < error_rate:
            # Random possibilities of deleting a base, changing a base to N, changing a base
            # to a random base, and inserting an N
            options = [lambda f,i: f[:i] + f[i+1:], lambda f,i: f[:i] + 'N' + f[i+1:], lambda f,i: f[:i] + random.choice("ATCG") + f[i+1:], lambda f,i: f[:i] + 'N' + f[i:]]
            fragment = random.choice(options)(fragment, random.randrange(0, len(fragment)))
        frags.append(fragment)
        if dna is seq:
            rs += length
        for i in range(rs, rs+rlen):
            coverage[i] += 1
        reads += 1
    return frags, reads, coverage

def matches(a, b, threshold):
    length = min(len(a), len(b))
    count = 0
    for i in range(length):
        if a[i] == b[i]:
            count += 1
    return count / length > threshold

def find_contigs(frags, threshold=0.75):
    contigs = set([])
    flen = len(frags)
    for i in range(flen):
        for j in range(flen):
            f1 = frags[i]
            f2 = frags[j]
            f1l = len(f1)
            f2l = len(f2)
            fminl = min(f1l, f2l)
            overlap = 0
            k = fminl - 1
            while k >= 1 and overlap == 0:
                if matches(f1[f1l-k:f1l], f2[0:k], threshold):
                    contig = f1[0:f1l-k] + f2
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                elif matches(f2[f2l-k:f2l], f1[0:k], threshold):
                    contig = f2[0:f2l-k] + f1
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                elif matches(complement(f1)[f1l-k:f1l], f2[0:k], threshold):
                    contig = complement(f1)[0:f1l-k] + f2
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                elif matches(complement(f2)[f2l-k:f2l], f1[0:k], threshold):
                    contig = complement(f2)[0:f2l-k] + f1
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                elif matches(f1[f1l-k:f1l], complement(f2)[0:k], threshold):
                    contig = f1[0:f1l-k] + complement(f2)
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                elif matches(f2[f2l-k:f2l], complement(f1)[0:k], threshold):
                    contig = f2[0:f2l-k] + complement(f1)
                    overlap = k
                    contigs.add((f1, f2, contig, overlap))
                k -= 1
    return contigs

seq = input('Sequence: ')
if len(seq) > 1 and seq[0] == '>':
    name = seq.rstrip()[1:].upper()
if not seq or seq[0] == '>':
    seq = ''
input_temp = 'ignored'
while input_temp:
    input_temp = input()
    if input_temp and input_temp[0] != '>':
        seq += input_temp.upper().rstrip()

fmin = int(input('Minimum fragment size [10]: ').strip() or 10)
fmax = int(input('Maximum fragment size [20]: ').strip() or 20)
fold = int(input('Coverage fold [30]: ').strip() or 30)
threshold = float(input('Overlap threshold [0.75]: ').strip() or 0.75)

frags, reads, coverage = build_fragments(seq, fmin, fmax, fold)

print("Number of sequence reads:", reads)
print("Coverage:", coverage)
print("Number of fragments containing N:", len([f for f in frags if 'N' in f]))

contigs = find_contigs(frags, threshold=threshold)

for contig in contigs:
    print('Frag1: {}, Frag2: {}, Contig: {}, Overlap: {}'.format(*contig))

print("{} total contigs.".format(len(contigs)))