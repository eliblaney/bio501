import sys
import re
import math

# Handy function to print a matrix nicely
# I threw this together earlier for the purpose of quick debugging, but I
# can also expand it out to a full function and write pretty code/comments if
# you would prefer that instead. It basically iterates each row of the matrix,
# transforms each row element into a string and rounds it to the nearest integer,
# pads it to 4 total characters, and joins all the rows together with nice borders.
print_matrix=lambda m:print('-'*(5*len(m[0])+3),*['| '+' '.join(map(lambda x:' '*(4-len(x))+x,map(lambda x:'{:.0f}'.format(x),r)))+' |'for r in m],'-'*(5*len(m[0])+3),sep='\n')if m else print('(empty matrix)')

def build_matrix(seq1, seq2, substitution_matrix, gap_score=-1, mismatch_score=0, match_score=1, align_type='global'):
    """Use the Needleman-Wunsch algorithm to build an alignment matrix for two sequences.

    Keyword arguments:
    seq1 -- A tuple of the first sequence name and its sequence.
    seq2 -- A tuple of the second sequence name and its sequence.
    substitution_matrix -- The substitution matrix to identify optimal replacements
    gap_score -- The penalty to apply to gaps in aligned sequences. Default is -1.
    mismatch_score -- The penalty to apply to gaps in aligned sequences. Default is 0, ignored for proteins.
    match_score -- The reward to apply to gaps in aligned sequences. Default is 1, ignored for proteins.
    align_type -- The type of alignment. Can be 'global', 'semiglobal', or 'local'.
    """
    # Quick lambdas to add corrections for global, semiglobal, and local alignment differences.
    # Global alignments will find the maximum score from the compared values in each cell,
    # semiglobal and local alignments will do the same but also have minimum values of zero
    # in the terminal gap regions. Local alignments will furthermore have a minimum value of
    # zero for all cells in the matrix.
    correct = lambda *scores: max(0, *scores) if align_type == 'local' else max(scores)
    terminal_correct = lambda *scores: max(0, *scores) if align_type != 'global' else max(scores)

    name1, seq1 = seq1
    name2, seq2 = seq2

    if len(seq1) < len(seq2):
        temp = (name1, seq1)
        name1, seq1 = (name1, seq1)
        name2, seq2 = temp

    N = len(seq1)
    M = len(seq2)

    matrix = []
    # Top-left corner zero to start the show
    matrix.append([0])
    # Top horizontal row (seq1)
    for i in range(1, N + 1):
        # Add the gap score to each successive cell
        val = matrix[0][i - 1] + gap_score
        # Correct it with zeros for semiglobal and local alignments
        matrix[0].append(terminal_correct(val))

    # Iterate rows (increasing seq2)
    for j in range(1, M + 1):
        # Side vertical row (seq2)
        # Add the gap score to each successive cell
        val = matrix[j - 1][0] + gap_score
        # Correct it with zeros for semiglobal and local alignments
        row = [terminal_correct(val)]
        seq2base = seq2[j - 1]
        # Iterate columns (increasing seq1)
        for i in range(1, N + 1):
            score = 0
            seq1base = seq1[i - 1]
            # Pseudocount for invalid amino acids like "X"
            score1 = matrix[j - 1][i - 1] + 1
            if substitution_matrix:
                if seq1base in substitution_matrix and seq2base in substitution_matrix[seq1base]:
                    score1 = matrix[j - 1][i - 1] + substitution_matrix[seq1base][seq2base]
            else:
                if seq1[i - 1] == seq2[j - 1]:
                    # Get the score for matched bases diagonally
                    score1 = matrix[j - 1][i - 1] + match_score
                else:
                    # Get the score for mismatched bases diagonally
                    score1 = matrix[j - 1][i - 1] + mismatch_score
            # Find the maximum value for substitution matrix values with
            # values coming from gaps in both directions, while also
            # correcting for local alignment minimum values.
            score = correct(
                score1,
                matrix[j - 1][i] + gap_score,
                row[i - 1] + gap_score
            )
            row.append(score)

        matrix.append(row)

    # Build the alignment state dictionary that keeps track of all the
    # user-provided parameters and current state of the alignments.
    return {
            'matrix': matrix,
            'substitution_matrix': substitution_matrix,
            'seq1': seq1,
            'seq2': seq2,
            'name1': name1,
            'name2': name2,
            'gap_score': gap_score,
            'mismatch_score': mismatch_score,
            'match_score': match_score,
            'align_type': align_type,
            'alignments': []
            }

def jukes_cantor(seq1, seq2, count_gaps=False):
    """Calculate the Jukes-Cantor distance for two sequences

    Keyword arguments:
    seq1 -- The first aligned sequence
    seq2 -- The second aligned sequence
    count_gaps -- Whether to count gaps as mismatches. Defaults to False.
    """
    difCtr = 0
    lenCtr = 0
    for k in range(len(seq1)):
        if seq1[k] != '-' and seq2[k] != '-':
            lenCtr += 1
            if seq1[k] != seq2[k]:
                difCtr += 1
        elif count_gaps and (seq1[k] == '-' or seq2[k] == '-'):
            lenCtr += 1
            difCtr += 1

    difCtr /= lenCtr
    return -3/4 * math.log(1 - 4/3 * difCtr)

def kimura(seq1, seq2):
    """Calculate the Kimuradistance for two sequences

    Keyword arguments:
    seq1 -- The first aligned sequence
    seq2 -- The second aligned sequence
    """
    transitions = 0
    transversions = 0
    total_substitutions = 0
    for i in range(len(seq1)):
        j = seq1[i]
        k = seq2[i]
        if j != '-' and k !='-':
            total_substitutions += 1
            if j != k:
                if j == 'A' and k == 'G' or j == 'G' and k == 'A':
                    transitions += 1
                elif j == 'T' and k == 'C' or j == 'C' and k == 'T':
                    transitions += 1
                else:
                    transversions += 1

    p = transitions / total_substitutions
    q = transversions / total_substitutions
    return -0.5 * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))
            
def scores(alignment, align1, align2, jukes_count_gaps=False):
    """Determine the identity score, gap score, and overall score for an alignment.

    Keyword arguments:
    alignment -- The alignment state dictionary, containing the gap penalties
    align1 -- The first ALIGNED sequence
    align2 -- The second ALIGNED sequence
    jukes_count_gaps -- Count gaps as mismatches for the Jukes-Cantor algorithm. Defaults to False.
    """
    gap_score = alignment['gap_score']
    
    identity = 0
    # Gaps are cumulatively added from both sequences
    gaps = align1.count('-') + align2.count('-')
    # The overall score is just the bottom right value in the matrix
    score = alignment['matrix'][-1][-1]
    # Calculate the number of matches
    for i in range(len(align1)):
        a = align1[i]
        b = align2[i]
        if a == b:
            identity += 1

    return identity, gaps, score, jukes_cantor(align1, align2, count_gaps=jukes_count_gaps), kimura(align1, align2)

def align_sequences(alignment, find_all=True, align1='', align2='', cur_row=None, cur_col=None, jukes_count_gaps=False):
    """Align some or all of the sequences provided from an alignment matrix.

    Keyword arguments:
    alignment -- The alignment state dictionary, containing the gap penalties
    find_all -- Whether to look for ALL optimal alignments or just one.
    align1 -- For internal use. The first ALIGNED sequence
    align2 -- For internal use. The second ALIGNED sequence
    cur_row -- For internal use. The current row of the matrix to examine.
    cur_col -- For internal use. The current column of the matrix to examine.
    """
    # Unpack the alignment state dictionary to variables.
    align_type, alignments, gap_score, match_score, matrix, mismatch_score, name1, name2, seq1, seq2, substitution_matrix = [v[1] for v in sorted(alignment.items())]
    alignments = []

    # The first time this function is called, we need to calculate where to start
    # finding the alignment.
    if cur_row == None or cur_col == None:
        cur_row = len(matrix) - 1
        cur_col = len(matrix[0]) - 1
        if align_type == 'semiglobal':
            max_val = None
            cur_row = len(matrix) - 1
            row = matrix[cur_row]
            num_cols = len(row)
            for j in range(num_cols):
                val = row[j]
            if max_val is None or val >= max_val:
                max_val = val
            cur_col = j
            num_gaps = num_cols - cur_col - 1
            align2 = '-' * num_gaps
            align1 = seq1[-num_gaps:]
        elif align_type == 'local':
            max_val = None
            for i in range(len(matrix)):
                row = matrix[i]
                for j in range(len(row)):
                    val = matrix[i][j]
                    if max_val is None or val >= max_val:
                        max_val = val
                        cur_row = i
                        cur_col = j

    # For semiglobal and global alignment, we can stop our recursion at the top-left
    # corner of the matrix. For local alignment, we stop when we find a zero.
    if cur_row <= 0 and cur_col <= 0 or (align_type == 'local' and matrix[cur_row][cur_col] == 0):
        # Reverse our aligned sequences and combine them with the calculated
        # scores, and then propagate back through the call stack.
        align1 = align1[::-1]
        align2 = align2[::-1]
        identity, gaps, score, jukes_cantor, kimura = scores(alignment, align1, align2, jukes_count_gaps)
        # Combine our alignment state dictionary with our new alignments, keeping
        # the old keys/values but overwriting the alignments item.
        return dict(alignment, alignments=[(align1, align2, identity, gaps, score, jukes_cantor, kimura)])

    # Make handy aliases for the different directions and cells that
    # will be referenced in determining potential paths to follow.
    up = matrix[cur_row - 1][cur_col]
    left = matrix[cur_row][cur_col - 1]
    diag = matrix[cur_row - 1][cur_col - 1]
    cur = matrix[cur_row][cur_col]

    # Helpful lambdas to make the following code more readable.
    # build will call the recursive function, copying over the
    # current state and retrieving only the alignments of the
    # subset.
    build = lambda a1, a2, ro, co: align_sequences(alignment, find_all, a1, a2, cur_row + ro, cur_col + co, jukes_count_gaps)['alignments']
    # Extend the alignments of this subsetted arrangement with
    # the contents of the recursive call, using a lambda so that
    # the following code is shorter and more readable.
    find_paths = lambda a, b, ro, co: alignments.extend(build(align1 + a, align2 + b, ro, co))

    # Try every possible direction to travel to build the sequences
    # We can keep track if we found a valid result to exit the function
    # early if we don't want to find all the alignments.
    found = False
    # Check horizontal gap movement
    if cur_row == 0 or left + gap_score == cur:
        seq1_base = seq1[cur_col - 1]
        # Recursively call the function with a smaller matrix
        find_paths(seq1_base, '-', 0, -1)
        found = True

    # Check vertical gap movement
    if (not found or find_all) and (cur_col == 0 or up + gap_score == cur):
        seq2_base = seq2[cur_row - 1]
        find_paths('-', seq2_base, -1, 0)
        found = True

    # Diagonal movement using substitution matrix
    if (not found or find_all) and (cur_row != 0 and cur_col != 0):
        seq1_base = seq1[cur_col - 1]
        seq2_base = seq2[cur_row - 1]
        if substitution_matrix:
            change = 1
            # Check substitution matrix diagonal movements
            if seq1_base in substitution_matrix and seq2_base in substitution_matrix[seq1_base]:
                change = substitution_matrix[seq1_base][seq2_base] 
            if diag + change == cur:
                find_paths(seq1_base, seq2_base, -1, -1)
                found = True
        else:
            matching = seq1_base == seq2_base
            if matching and diag + match_score == cur:
                find_paths(seq1_base, seq2_base, -1, -1)
                found = True
            # Check mismatching movements
            elif not matching and diag + mismatch_score == cur:
                find_paths(seq1_base, seq2_base, -1, -1)
                found = True

    # Combine the alignments into the state dictionary
    return dict(alignment, alignments=alignments)

def print_alignment(alignment, substitution_matrix_name=None, num=1, output_buffer=print, print_jukes=True, print_kimura=True):
    """Output the alignment for paired sequences to a buffer.

    Keyword arguments:
    alignment -- The alignment state matrix.
    num -- The number of sequence pairs to print alignments for. Defaults to 1.
    output_buffer -- The buffer to send the output to. Defaults to the standard print function.
    print_jukes -- Whether to print the Jukes-Cantor distance. Defaults to True.
    print_kimura -- Whether to print the Kimura distance. Defaults to True.
    """
    # If the output buffer isn't print, then we need to manually
    # add a newline. Print will handle it automatically, so no
    # extra modifications are needed for the default argument.
    if output_buffer != print:
        old_buffer = output_buffer
        # Call the same buffer that was passed in, but add a newline
        output_buffer = lambda string: old_buffer(string + '\n')

    # Grab the aligned sequences from the state matrix
    alignments = alignment['alignments']

    # If a number argument wasn't
    if num == 0:
        num = len(alignments)

    # If we have too many alignments than ones that we
    # want to print, we want to reduce what is present.
    while len(alignments) > num:
        # Find the alignment pair with the worst score,
        # and remove each one until we have an acceptable
        # number of alignment pairs.
        # Because all the alignments in this algorithm
        # have the same score, this sorting does not have
        # any particular effect. However, if the scoring
        # metric changes in some way, then this is a
        # useful inclusion to the program.
        worst_score = None
        worst_index = len(alignments)
        for i in range(len(alignments)):
            a = alignments[i]
            if worst_score is None or a[4] < worst_score:
                worst_score = a[4]
                worst_index = i
        alignments = alignments[:i] + alignments[i+1:]

    # Print out all the alignments and their various scores.
    for align1, align2, identity, gaps, score, jukes_cantor, kimura in alignments:
        # Determine the string length of the characters in the
        # sequence, so that the indices can be aligned at the
        # beginning of each aligned segment for each line.
        length = len(align1)
        max_prefix_len = len(str(length)) + 1
        
        # Print out the header with formatted parameters
        output_buffer('#===================================')
        output_buffer("# Aligned sequences: 2")
        output_buffer('# Parameters:')
        if substitution_matrix_name:
            output_buffer('#   - sequence type: protein')
            output_buffer('#   - substitution matrix: {}'.format(substitution_matrix_name))
        else:
            output_buffer('#   - sequence type: nucleotide')
            output_buffer('#   - match score: {}'.format(alignment['match_score']))
            output_buffer('#   - mismatch score: {}'.format(alignment['mismatch_score']))
        output_buffer('#   - gap score: {}'.format(alignment['gap_score']))
        output_buffer('#   - alignment type: {}'.format(alignment['align_type']))
        output_buffer('#')
        output_buffer('# Length: {}'.format(length))
        output_buffer('# Identity: {}/{} ({:.1f}%)'.format(identity, length, identity/length*100))
        output_buffer('# Gaps: {}/{} ({:.1f}%)'.format(gaps, length, gaps/length*100))
        output_buffer('# Score: {:.1f}'.format(score))
        if print_jukes:
            output_buffer('# Jukes-Cantor: {:.4f}'.format(jukes_cantor))
        if print_kimura:
            output_buffer('# Kimura: {:.4f}'.format(kimura))
        output_buffer('#===================================')
        output_buffer('')

        line_wrap = 60
        # Wrap the aligned sequences along 60 character intervals
        for i in range(0, length, line_wrap):
            # The prefix is just the sequence position number
            prefix = str(i + 1)
            # The prefix length
            pl = len(prefix)
            # Pad the prefix with spaces to get to a standard length
            l1 = prefix + ' ' * (max_prefix_len - pl)
            # Skip the line number, just pad to the max length
            l2 = ' ' * max_prefix_len
            # Copy the first line's prefix
            l3 = l1
            # The final sequence position can hit the end
            # of the sequence
            end = min(i + line_wrap, length)

            # Iterate through the truncated line to print out the
            # alignments between each sequence
            for j in range(i, end):
                a = align1[j]
                b = align2[j]
                l1 += a
                l3 += b
                if a == b:
                    l2 += '|'
                elif a == '-' or b == '-':
                    l2 += ' '
                else:
                    l2 += '.'
            output_buffer(l1 + ' ' + str(end))
            output_buffer(l2)
            output_buffer(l3 + ' ' + str(end))
            output_buffer('')

def align_multiple(seqs, smatrix, gap_score=-1, mismatch_score=0, match_score=1, align_type='global'):
    """Align multiple sequences.

    Keyword arguments:
    seqs -- An array of sequences to align
    gap_score -- The penalty to apply to gaps in aligned sequences. Default is -1.
    mismatch_score -- The penalty to apply to mismatches in aligned sequences. Default is 0, ignored for proteins.
    match_score -- The penalty to apply to matches in aligned sequences. Default is 1, ignored for proteins.
    align_type -- The type of alignment. Can be 'global', 'semiglobal', or 'local'.
    """
    # Get optimal alignments for every combination of sequences
    # Intensive process, but it will allow us to select the best reference sequence.
    aggregated_alignments = [[align_sequences(build_matrix(refseq, seq, smatrix, gap_score, align_type=align_type), find_all=False) for seq in seqs if seq != refseq] for refseq in seqs]

    # Find the best alignment state matrix having the lowest number of mismatches
    best_alignments = None
    best_mismatches = None
    for alignments in aggregated_alignments:
        mismatches = 0
        for alignment in alignments:
            refseq, subseq, _, _, _ = alignment['alignments'][0]
            # Count mismatches, including gaps as mismatches
            for i in range(len(subseq)):
                if subseq[i] != refseq[i]:
                    mismatches += 1
        if best_mismatches == None or mismatches < best_mismatches:
            best_mismatches = mismatches
            best_alignments = alignments

    # This best alignment may not work perfectly for all sequences,
    # so we will quickly pad any weird shortened sequences using gaps.
    for alignment in best_alignments:
        refseq, subseq, _, _, _ = alignment['alignments'][0]
        lrefseq = len(refseq)
        lsubseq = len(subseq)
        if lrefseq > lsubseq:
            alignment['alignments'][0]['align2'] = subseq + ('-' * (lrefseq - lsubseq))

    # We now have the best alignments: all are being compared with one reference sequence,
    # and that combination has the lowest overall number of mismatches/gaps.
    return best_alignments

# Adapted from previous homework assignment
def short_name(name, max_chars=16):
    """Creates a shortened name based on the header or sequence number.

    Keyword arguments:
    name -- The long name of the sequence or the sequence number.
    max_chars -- The number of characters to constrain the name to.
                 Default to 16. If 0 or below, ignored.
    """
    # If the name is actually an integer, that means there was no
    # FASTA header, and the name is just a sequence number. We will
    # return something more meaningful by prepending "Sequence "
    if isinstance(name, int):
        return "Sequence {}".format(name)
    # If our name has exceeded the character limit, we need to
    # return a string that splices out the middle characters of
    # our name and replaces it with an ellipsis so that the length
    # of our name is exactly equal to the character limit.
    if max_chars > 0 and len(name) > max_chars:
        return name[:max_chars - 5] + "..." + name[-2:]
    # Otherwise, the name is shorter than the character limit,
    # so it is completely valid to return on its own.
    return name

def print_multiple_alignments(alignments, substitution_matrix_name=None, output_buffer=print, line_wrap=60, prefix_len=16, print_jukes=True, print_kimura=True):
    """Output the multiple alignment for a series of sequences.

    Keyword arguments:
    alignments -- A list of alignment state matrices.
    output_buffer -- The buffer to send the output to. Defaults to the standard print function.
    line_wrap -- The number of characters on each line
    prefix_len -- The maximum number of characters at the beginning of each line in the alignment
    print_jukes -- Whether to print the Jukes-Cantor distance. Defaults to True.
    print_kimura -- Whether to print the Kimura distance. Defaults to True.
    """
    # If the output buffer isn't print, then we need to manually
    # add a newline. Print will handle it automatically, so no
    # extra modifications are needed for the default argument.
    if output_buffer != print:
        old_buffer = output_buffer
        # Call the same buffer that was passed in, but add a newline
        output_buffer = lambda string: old_buffer(string + '\n')

    # Use the first_alignment as a shortcut for the scores and refseq
    first_alignment = alignments[0]
    # Length corresponds to the length of the aligned refseq
    length = len(first_alignment['alignments'][0][0])
    # Print out the header with formatted parameters
    output_buffer('#===================================')
    output_buffer("# Aligned sequences: {}".format(len(alignments) + 1))
    output_buffer('# Parameters:')
    if substitution_matrix_name:
        output_buffer('#   - sequence type: protein')
        output_buffer('#   - substitution matrix: {}'.format(substitution_matrix_name))
    else:
        output_buffer('#   - sequence type: nucleotide')
        output_buffer('#   - match score: {}'.format(alignment['match_score']))
        output_buffer('#   - mismatch score: {}'.format(alignment['mismatch_score']))
    output_buffer('#   - gap score: {}'.format(first_alignment['gap_score']))
    output_buffer('#   - alignment type: {}'.format(first_alignment['align_type']))
    output_buffer('#')
    output_buffer('# RefSeq: {}'.format(first_alignment['name1']))
    output_buffer('# Length: {}'.format(length))
    # Print out the scores for each of the aligned sequences
    i = 2
    for alignment in alignments:
        _, _, identity, gaps, score, jukes_cantor, kimura = alignment['alignments'][0]
        name = short_name(alignment['name2'], 10)
        output_buffer('# [{}] Identity: {}/{} ({:.1f}%)'.format(name, identity, length, identity/length*100))
        output_buffer('# [{}] Gaps: {}/{} ({:.1f}%)'.format(name, gaps, length, gaps/length*100))
        output_buffer('# [{}] Score: {:.1f}'.format(name, score))
        if print_jukes:
            output_buffer('# [{}] Jukes-Cantor: {:.4f}'.format(jukes_cantor))
        if print_kimura:
            output_buffer('# [{}] Kimura: {:.4f}'.format(kimura))
        i += 1
    output_buffer('#===================================')
    output_buffer('')
    
    # A little different than the other method of printing alignments,
    # this code will build partitions for each of the sequences so that they
    # can all be printed in chunks simultaneously.
    partitions = []
    num_lines = 0
    for alignment in alignments:
        refseq, aligned_seq, identity, gaps, score, jukes_cantor, kimura = alignment['alignments'][0]
        # On the first iteration, add the refseq so that it is first in the list
        if not partitions:
            # Take advantage of regex to make an array of subsequences up to 60 chars
            refseq_partition = re.findall(".{1," + str(line_wrap) + "}", refseq)
            # All sequences will necessarily have the same number of lines, so we
            # only need to do this calculation once.
            num_lines = len(refseq_partition)
            # Get the shortened name and store it with the partitions as a tuple
            name = short_name(alignment['name1'], prefix_len-2)
            partitions.append((name, refseq_partition))
        # Same as above, but with each sequence
        name = short_name(alignment['name2'], prefix_len-2)
        partitions.append((name, re.findall(".{1," + str(line_wrap) + "}", aligned_seq)))

    # For every line, we can print the subsequences that correspond with that line.
    # We also figure out what is matching, so we can print appropriate symbols.
    for line in range(num_lines):
        mismatches = []
        refseq = partitions[0][1][line]
        for name, partition in partitions:
            subseq = partition[line]
            # Find mismatches among ALL sequences and keep track of their indices
            for i in range(len(subseq)):
                if i > len(refseq) or (not i in mismatches and subseq[i] != refseq[i]):
                    mismatches.append(i)
            padding = prefix_len - len(name) - 2
            # Print shortened name, with extra padding up to prefix_len
            output_buffer(name + ': ' + (' ' * padding) + subseq)

        # By default, everything is matching with "*"
        special = ' ' * prefix_len + '*' * len(refseq)
        # Then we replace the indices of mismatches with spaces.
        for i in range(len(mismatches)):
            index = mismatches[i] + prefix_len
            # Replace with space at the mismatch index
            special = special[:index] + ' ' + special[index+1:]
        # Output the symbols and go on to the next line
        output_buffer(special)
        output_buffer('')

# Get AMINO ACID sequences from user input
seqs = []
while True:
    name = len(seqs) + 1
    seq = input('Sequence {}: '.format(name))
    if len(seq) > 1 and seq[0] == '>':
        name = seq.rstrip()[1:].upper()
    if not seq or seq[0] == '>':
        seq = ''
    input_temp = 'ignored'
    while input_temp:
        input_temp = input()
        if input_temp and input_temp[0] != '>':
            seq += input_temp.upper().rstrip()

    if seq:
        seqs.append((name, seq))
    else:
        break

# Get parameters from user input
smatrix = input('Substitution Matrix [nt]: ').lower()
if smatrix:
    smatrix += ".txt"
else:
    smatrix = None
gap_score = input('Gap score [-1]: ')
if gap_score:
    gap_score = float(gap_score)
else:
    gap_score = -1

mismatch_score = 0
match_score = 1
if smatrix == None:
    mismatch_score_str = input('Mismatch score [0]: ')
    if mismatch_score_str:
        mismatch_score = float(mismatch_score_str)
    match_score_str = input('Match score [1]: ')
    if match_score_str:
        match_score = float(match_score_str)

align_type = input('Align type [global]: ')
if align_type:
    align_type = align_type.strip().lower()
else:
    align_type = 'global'

print_jukes_cantor_kimura = input('Stats (1 = Jukes-Cantor, 2 = Kimura, 3 = both) [3]: ')
if print_jukes_cantor_kimura:
    print_jukes_cantor_kimura = int(print_jukes_cantor_kimura.strip())
else:
    print_jukes_cantor_kimura = 3

amino_acids = None
smatrix_name = None
if smatrix:
    smatrix_file = open(smatrix)
    smatrix = [row.split(",") for row in smatrix_file.read().split("\n") if row]
    smatrix_name = smatrix[0][0][1:]
    smatrix_file.close()

    amino_acids = smatrix[1]
    smatrix = [dict(zip(amino_acids, map(float, row))) for row in smatrix[2:]]
    smatrix = {acid: smatrix[i] for i, acid in enumerate(amino_acids)}

# I used recursion, but Python's recursion limit is quite
# low. Unfortunate, but this is something we can change. Here,
# I dynamically change the recursion limit based on the length
# of our sequences so that we don't run into any problems.
sys.setrecursionlimit(max(1000, *map(len, seqs)) * 10)

num_seqs = len(seqs)
if num_seqs <= 1:
    print("Not enough sequences to align.")

alignments = []
if num_seqs == 2:
    do_print_matrix = 'y' == input('Print matrix [n]? ').lower()
    num_alignments = input('Number of alignments to show [1]: ')
    if num_alignments and num_alignments.lower() == 'all':
        num_alignments = 0
    elif num_alignments:
        num_alignments = int(num_alignments)
    else:
        num_alignments = 1

    # Do the work by passing in our arguments!
    # Create a new alignment state dictionary
    alignment = build_matrix(seqs[0], seqs[1], smatrix, gap_score, mismatch_score, match_score, align_type=align_type)
    # Align the sequences quite easily just by passing in
    # the alignment state dictionary. We can choose to reduce
    # the information that we get by passing in False to find_all,
    # limiting the number of optimal sequences obtained to just one.
    alignment = align_sequences(alignment, num_alignments != 1)
    if do_print_matrix:
        print_matrix(alignment['matrix'])

    # Print the alignment easily by passing in the alignment
    # state dictionary.
    print_alignment(alignment, substitution_matrix_name=smatrix_name, num=num_alignments, print_jukes=print_jukes_cantor_kimura % 2 == 1, print_kimura=print_jukes_cantor_kimura >= 2)
    alignments = [alignment]
else:
    alignments = align_multiple(seqs, smatrix, gap_score, mismatch_score, match_score, align_type=align_type)
    print_multiple_alignments(alignments, substitution_matrix_name=smatrix_name, print_jukes=print_jukes_cantor_kimura % 2 == 1, print_kimura=print_jukes_cantor_kimura >= 2)