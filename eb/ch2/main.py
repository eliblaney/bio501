### Eli Blaney ###
###  BIO 501   ###
### Chapter 2  ###
### 2022-02-10 ###

def short_name(names, num, max_chars=10):
    """Creates a shortened name based on the header or sequence number.

    Keyword arguments:
    names -- The list of names which contains all names for the sequences
    num -- The sequence number corresponding with the name. NOT the index.
    max_chars -- The number of characters to constrain the name to.
                 Default to 10. If 0 or below, ignored.
    """
    if not names:
        return "Seq #{}".format(num)
    # Get the name associated with our sequence. The sequence number
    # is just the index plus one, so we can decrement to get it.
    name = names[num - 1]
    # If the name is actually an integer, that means there was no
    # FASTA header, and the name is just a sequence number. We will
    # return something more meaningful by prepending "Seq #"
    if isinstance(name, int):
        return "Seq #{}".format(name)
    # If our name has exceeded the character limit, we need to
    # return a string that splices out the middle characters of
    # our name and replaces it with an ellipsis so that the length
    # of our name is exactly equal to the character limit.
    if max_chars > 0 and len(name) > max_chars:
        return name[:max_chars - 5] + "..." + name[-2:]
    # Otherwise, the name is shorter than the character limit,
    # so it is completely valid to return on its own.
    return name

def transcribe(dnas):
    """Transcribe a list of DNA sequences to RNA. Maintains orientation.

    Keyword arguments:
    dnas -- A list of DNA sequences.
    """
    rnas = []
    for seq in dnas:
        # To transcribe, just replace thymine with uracil
        rnas.append(seq.replace('T', 'U'))

    return rnas

def translate(rnas, offset=0):
    """Translate a list of RNA sequences to amino acids.

    Keyword arguments:
    rnas -- A list of mRNA sequences from the coding DNA
            strand, in the 5' to 3' orientation.
    offset -- An optional reading frame adjustment. Default is 0.
    """
    # Useful codon table from a previous assignment
    codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'U', 'ACC':'U', 'ACG':'U', 'ACU':'U',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'X', 'UAG':'X',
        'UGC':'C', 'UGU':'C', 'UGA':'X', 'UGG':'W'
    }

    proteins = []
    for seq in rnas:
        aa_seq = ""
        # We iterate through the sequence starting from
        # the beginning of the mRNA, adjusted by the
        # offset parameter. We jump 3 bases each time so
        # that we can select entire codons without doing
        # too much calculation.
        for i in range(3 + offset, len(seq), 3):
            # Select the codon by finding the last 3
            # bases in our sequence that we've reached.
            codon = seq[i-3:i]
            # If there are extra bases at the end, we don't
            # want them. Just translate full sets of 3 bases.
            if len(codon) == 3:
                # Grab the correct amino acid and add it
                # to our amino acid string sequence.
                aa_seq += codon_table[codon]

        # Add each amino acid sequence to our proteins list
        proteins.append(aa_seq)
    return proteins

def align(seq1, seq2, warn_frameshifts=False):
    """Align a set of sequences to a given reference sequence.
    Constrained to textbook assumptions that the sequences have
    only had zero or one mutational events, causing consecutive
    gaps in only one of the sequences.

    Keyword arguments:
    seq1 -- The first sequence to align.
    seq2 -- The second sequence to align.
    warn_frameshifts -- Whether to identify nucleotide sequences that
                        contain a frameshift mutation. Default to False.
    """
    # We don't need to align anything if they're the same
    # length, so return early in that case.
    seq1_length = len(seq1)
    seq2_length = len(seq2)
    if seq1_length == seq2_length:
        return seq1, seq2

    # From this point forward, we want to assume that seq1
    # is the longer strand. Therefore, if seq2 is longer,
    # we will swap our variables.
    # args_reversed will keep track of this switch so that
    # the returned values can maintain the order that they
    # were given in.
    args_reversed = False
    if seq1_length < seq2_length:
        # Reusing args_reversed as a temporary storage
        # variable; it will evaluate truthy, which is what
        # we want to happen. Python's dynamic typing allows
        # this to work well, saving a minor amount of memory.
        args_reversed = seq1_length
        seq1_length = seq2_length
        seq2_length = args_reversed
        args_reversed = seq2
        seq2 = seq1
        seq1 = args_reversed

    # Primitive algorithm which the textbook wants to use to
    # determine the indel gaps; subtract their lengths to find
    # their number and compare each possibility of indel
    # combination to find the best option. Inefficient due
    # to polynomial runtime, ignores mixed insertions and
    # deletions.
    num_gaps = seq1_length - seq2_length

    # If we enable the warn_frameshifts flag, we know that
    # we are dealing with nucleotides, as indicated by the
    # calling scope. Because of this, we can indicate that
    # there is a frameshift mutation with certainty if there
    # is a gap number that is not a multiple of three.
    # We indicate this frameshift by returning None for the
    # expected aligned sequence. Ideally, we would like to
    # properly handle this frameshift, but the assignment
    # just wants us to identify them until we cover better
    # tools to do so.
    if warn_frameshifts and num_gaps % 3 != 0:
        if args_reversed:
            return None, seq1
        return seq1, None

    # Now we want to brute-force all possible combinations
    # of consecutive alignments. We can do this by iterating
    # through every index and inserting consecutive dashes,
    # comparing them to find the best alignment with the
    # least number of mutations.
    best_seq2 = None
    best_seq2_mutations = -1
    # Shaving off a couple milliseconds from this algorithm by
    # doing this string multiplication only once and reusing
    # the variable.
    dashes = "-" * num_gaps
    for i in range(seq2_length):
        # Insert the consecutive dashes at each index
        modified_seq2 = seq2[:i] + dashes + seq2[i:]

        # Compare the reference sequence with this modified
        # sequence to find how many mutations are present.
        # We grab the first index because the function returns
        # a list of aggregated mutations, and we only pass in
        # two sequences to compare, so there is only one
        # mutation list to look at. do_align is set to False
        # because we want to avoid an infinite loop in some
        # edge cases.
        mutations = compare([seq1, modified_seq2], do_align=False)[0]
        # The compare function may return None if there is some
        # unexpected frameshift mutation, in which case we know
        # it is definitely not the best, so we skip it.
        if mutations == None:
            continue
        # Find the number of mutations so we can compare it with
        # the best aligned sequence that we've found so far.
        mutations = len(mutations)
        if best_seq2_mutations == -1 or mutations < best_seq2_mutations:
            best_seq2 = modified_seq2
            best_seq2_mutations = mutations

    # We kept track if we needed to reverse the function
    # arguments, so now we can return them in the same order
    # that the calling scope provided, so they can maintain
    # safety in assigning new variables.
    if args_reversed:
        return best_seq2, seq1
    return seq1, best_seq2

def compare(seqs, do_align=True, warn_frameshifts=False):
    """Compare a list of sequences to find all possible mismatch
       mutations. Returns an aggregated list of mutations.

    Keyword arguments:
    seqs -- The sequences to compare. Can be any sequence containing
            single-character elements, such as DNA, RNA, or amino acids.
            The first element in the list is considered the reference
            sequence; in other words, all other sequences will be
            compared against it.
    do_align -- Whether to align each sequence before comparing them. Default
                to True.
    warn_frameshifts -- Whether to identify nucleotide sequences that
                        contain a frameshift mutation. Default to False.
    """
    # We need at least a reference sequence and something
    # to compare with it.
    if len(seqs) < 2:
        return None

    # Our first sequence is the reference sequence.
    refseq = seqs[0]
    length = len(refseq)

    seq_num = 1
    # We can keep track of all the mutations in a 2D array,
    # where we aggregate each sequence's list of mutations.
    aggregated_mutations = []
    # We can go through all the sequences, skipping the very
    # first sequence which is the reference.
    for seq in seqs[1:]:
        seq_num += 1
        # We will avoid supporting sequence alignment for now,
        # until we implement it later.
        aligned_refseq = refseq
        if len(seq) != length:
            if do_align:
                aligned_refseq, seq = align(refseq, seq, warn_frameshifts)
                if seq == None:
                    # Skip the comparison for this element by
                    # "continuing" on to the next iteration in the
                    # for loop
                    aggregated_mutations.append(None)
                    continue
            else:
                # Same as above comment
                aggregated_mutations.append(None)
                continue

        mutations = []
        for i in range(length):
            # Iterate through all bases in both sequences
            ref_element = aligned_refseq[i]
            new_element = seq[i]
            # Compare each element at each index in both strings
            # If they are different, report that.
            if ref_element != new_element:
                # The mutation can be represented by
                # putting the reference sequence element, the
                # position of the mutation, and the mutated
                # element. Deletions can be represented with
                # a delta (∆).
                if new_element == "-":
                    mutation = "∆{}{}".format(ref_element, i + 1)
                else:
                    mutation = "{}{}{}".format(ref_element, i + 1, new_element)
                mutations.append(mutation)
        aggregated_mutations.append(mutations)

    return aggregated_mutations

def print_compare(seqs, names=None, do_align=True, warn_frameshifts=False):
    """Compare a list of sequences to find all possible mismatch
       mutations. Returns nothing; prints the found results.

    Keyword arguments:
    names -- The names of each sequence, for more helpful printing. Defaults
             to using sequence numbers in place of names.
    seqs -- The sequences to compare. Can be any sequence containing
            single-character elements, such as DNA, RNA, or amino acids.
            The first element in the list is considered the reference
            sequence; in other words, all other sequences will be
            compared against it.
    do_align -- Whether to align each sequence before comparing them. Defaults
                to True.
    warn_frameshifts -- Whether to identify nucleotide sequences that
                        contain a frameshift mutation. Defaults to False.
    """
    # We need at least a reference sequence and something
    # to compare with it.
    if len(seqs) < 2:
        print("Not enough sequences to compare.")
        return

    # Our first sequence is the reference sequence.
    refseq = seqs[0]
    print("Reference sequence:", short_name(names, 1, 0))
    aggregated_mutations = compare(seqs, align, warn_frameshifts)

    for i in range(1, len(seqs)):
        seq = seqs[i]
        mutations = aggregated_mutations[i - 1]
        name = short_name(names, i + 1)
        if mutations == None:
            if warn_frameshifts and do_align:
                print("[{}] Frameshift mutation present. (length = {}, refseq = {})".format(name, len(seq), len(refseq)))
            elif warn_frameshifts:
                print("[{}] Sequences unaligned or frameshift mutation present. (length = {}, refseq = {})".format(name, len(seq), len(refseq)))
            else:
                print("[{}] Sequences unaligned; incompatible lengths. (length = {}, refseq = {})".format(name, len(seq), len(refseq)))
        elif mutations:
            print("[{}] {} // {} total mutations found.".format(name, ', '.join(mutations), len(mutations)))
        else:
            print("[{}] no mismatches found - strings are identical".format(name))

# First, we need to find out what input files contain the 
# sequences that we are interested in. We'll accept
# them as paths from the user.
inputs = []
user_input = "ignored"
# Until the user pressees enter, we can continue to accept
# files, because "\n" is falsey, but file paths are truthy
while user_input:
    user_input = input("Enter input file paths (or enter to stop): ")
    if user_input:
        inputs.append(user_input)

# Determine from the user if we want to compare DNA, RNA,
# or amino acids. The numbering system works fine to avoid
# the issue of dealing with various ways that the user could
# represent their intentions.
compare_flag = input("Enter 1-4 to compare DNA, RNA, Amino Acids, or both DNA and Amino Acids [4]: ");
# We want this variable as an integer
if compare_flag:
    compare_flag = int(compare_flag)
else:
    # If they didn't give anything, or it's invalid,
    # we will default to 4 (DNA and amino acid comparisons)
    compare_flag = 4

# This names list will hold the FASTA header lines that
# correspond to each DNA sequence. If no FASTA header line
# beginning with ">" is found, the number of the sequence
# will be used in its place.
names = []
# This dnas list holds the DNA sequences that are found in
# each file, collected together. We may find RNA sequences,
# but they will be reverse-transcribed back to DNA later.
dnas = []
for filename in inputs:
    file = open(filename)
    seq = ""
    for line in file:
        # Make sure to remove any newline characters
        temp_seq = line.rstrip()
        # If we find a ">", we know that there is a
        # new sequence in the file.
        if temp_seq[0] == '>':
            # We should record the previous sequence
            # (if there is one) in our dnas list and
            # reset our seq variable for the new sequence.
            if seq:
                dnas.append(seq)
                seq = ""
            # If the first sequence in the file is unnamed
            # for some reason, but we encounter a named
            # sequence later, we should record the initial
            # unnamed sequence and label it with the sequence
            # number.
            if len(dnas) > len(names):
                names.append(len(dnas))
            # Finally, record the name of the new sequence
            # based on the header that we found.
            names.append(temp_seq[1:])
            # Because this line will not have a sequence
            # on the same line as the header, we can skip
            # the rest of the code in our loop and "continue"
            # on to the next line in the file.
            continue
        
        # Until we get to the end of the file or find a
        # new sequence, we just attach the new line of
        # bases to the existing string.
        seq += temp_seq.upper()

    # Record the final sequence in the file, if there is one
    if seq:
        dnas.append(seq)
        # Again, ensure it has a name. If not, give it the
        # default name of its sequence number.
        if len(dnas) > len(names):
            names.append(len(dnas))
    file.close()


# Now, we need to ensure that all the sequences are uniformly
# all DNA, oriented 5' to 3', and coding strands. We do this
# by asking the user and making necessary corrections.
for i in range(len(dnas)):
    seq = dnas[i]
    # Show the user the short name of the sequence so that they
    # know which sequence we are asking them about.
    print("-----------------------")
    print(":: {}".format(short_name(names, i + 1, -1)))

    # If the sequence is RNA, we need to transform uracil to thymine
    is_rna = "rna" == input("Enter DNA or RNA for sequence type [DNA]: ").lower()
    if is_rna:
        dnas[i] = seq.upper().replace('U', 'T')

    # If the sequence is 3' to 5', we need to reverse the sequence
    is_three_prime = "3-5" == input("Enter 3-5 or 5-3 for DNA orientation [5-3]: ")
    if is_three_prime:
        dnas[i] = seq[::-1]

    # If the sequence is a template strand, we find the complement to
    # get the coding strand, also in the 5' to 3' direction.
    is_template = "template" == input("Enter template or nontemplate for strand type [nontemplate]: ")
    if is_template:
        dnas[i] = seq[::-1].replace('T', 'a').replace('A' 't').replace('G', 'c').replace('C', 'g')

    # Now, all our sequences are corrected and ready for
    # transcription, translation, and comparisons.

# Transcribe and translate for our comparisons
rnas = transcribe(dnas)
amino_acids = translate(rnas)

# Figure out what we want to compare and display
# the correct data. In the case that compare_flag
# is equal to 4, we will print out both the
# DNA comparison and amino acid comparison.
if compare_flag == 1 or compare_flag == 4:
    print("-----------------------")
    print("     Comparing DNA     ")
    print("-----------------------")
    print_compare(dnas, names, do_align=True, warn_frameshifts=True)
if compare_flag == 2:
    print("-----------------------")
    print("     Comparing RNA     ")
    print("-----------------------")
    print_compare(rnas, names, do_align=True, warn_frameshifts=True)
if compare_flag >= 3:
    print("-----------------------")
    print(" Comparing Amino Acids ")
    print("-----------------------")
    print_compare(amino_acids, names, do_align=True, warn_frameshifts=False)
