import sys
sys.setrecursionlimit(10000)

def build_matrix(seq1, seq2, gap_score=-1, mismatch_score=0, match_score=1, align_type='global'):
    correct = lambda *scores: max(0, *scores) if align_type == 'local' else max(scores)
    terminal_correct = lambda *scores: max(0, *scores) if align_type != 'global' else max(scores)
    N = len(seq1)
    M = len(seq2)

    matrix = []
    matrix.append([0])
    # Top horizontal row (seq1)
    for i in range(1, N + 1):
        val = matrix[0][i - 1] + gap_score
        matrix[0].append(terminal_correct(val))

    # Iterate rows (increasing seq2)
    for j in range(1, M + 1):
        # Side vertical row (seq2)
        val = matrix[j - 1][0] + gap_score
        row = [terminal_correct(val)]
        # Iterate columns (increasing seq1)
        for i in range(1, N + 1):
            score = 0
            if seq1[i - 1] == seq2[j - 1]:
                score = matrix[j - 1][i - 1] + match_score
            else:
                score = matrix[j - 1][i - 1] + mismatch_score
            score = correct(score, matrix[j - 1][i] + gap_score, row[i - 1] + gap_score)
            row.append(score)

        matrix.append(row)

    return {
            'matrix': matrix,
            'seq1': seq1,
            'seq2': seq2,
            'gap_score': gap_score,
            'mismatch_score': mismatch_score,
            'match_score': match_score,
            'align_type': align_type,
            'alignments': []
            }

def scores(alignment, align1, align2):
    gap_score = alignment['gap_score']
    mismatch_score = alignment['mismatch_score']
    match_score = alignment['match_score']
    
    identity = 0
    gaps = align1.count('-') + align2.count('-')
    score = 0.0
    for i in range(len(align1)):
        a = align1[i]
        b = align2[i]
        if a == b:
            identity += 1
            score += match_score
        elif a == '-' or b == '-':
            score += gap_score
        else:
            score += mismatch_score

    return identity, gaps, score

def align_sequences(alignment, find_all=True, align1='', align2='', cur_row=None, cur_col=None):
    align_type, alignments, gap_score, match_score, matrix, mismatch_score, seq1, seq2 = [v[1] for v in sorted(alignment.items())]
    alignments = []

    if cur_row == None or cur_col == None:
        cur_row = len(matrix) - 1
        cur_col = len(matrix[0]) - 1
        if align_type == 'local':
            max_val = None
            for i in range(len(matrix)):
                row = matrix[i]
                for j in range(len(row)):
                    val = matrix[i][j]
                    if max_val is None or val > max_val:
                        max_val = val
                        cur_row = i
                        cur_col = j

    if cur_row == 0 and cur_col == 0 or (align_type == 'local' and matrix[cur_row][cur_col] == 0):
        align1 = align1[::-1]
        align2 = align2[::-1]
        identity, gaps, score = scores(alignment, align1, align2)
        return dict(alignment, alignments=[(align1, align2, identity, gaps, score)])

    up = matrix[cur_row - 1][cur_col]
    left = matrix[cur_row][cur_col - 1]
    diag = matrix[cur_row - 1][cur_col - 1]
    cur = matrix[cur_row][cur_col]
    seq1_base = seq1[cur_col - 1]
    seq2_base = seq2[cur_row - 1]
    matching = seq1_base == seq2_base

    build = lambda a1, a2, ro, co: align_sequences(alignment, find_all, a1, a2, cur_row + ro, cur_col + co)['alignments']
    find_paths = lambda a, b, ro, co: alignments.extend(build(align1 + a, align2 + b, ro, co))

    found = False
    if cur_row == 0 or left + gap_score == cur:
        find_paths(seq1_base, '-', 0, -1)
        found = True

    if (not found or find_all) and (cur_col == 0 or up + gap_score == cur):
        find_paths('-', seq2_base, -1, 0)
        found = True

    if (not found or find_all) and (cur_row != 0 and cur_col != 0):
        if matching and diag + match_score == cur:
            find_paths(seq1_base, seq2_base, -1, -1)
        elif not matching and diag + mismatch_score == cur:
            find_paths(seq1_base, seq2_base, -1, -1)

    return dict(alignment, alignments=alignments)

def print_alignment(alignment, num=1, output_buffer=print):
    if output_buffer != print:
        old_buffer = output_buffer
        output_buffer = lambda string: old_buffer(string + '\n')
    alignments = alignment['alignments']
    if num == 0:
        num = len(alignments)
    while len(alignments) > num:
        worst_score = None
        worst_index = len(alignments)
        for i in range(len(alignments)):
            a = alignments[i]
            if worst_score is None or a[4] < worst_score:
                worst_score = a[4]
                worst_index = i
        alignments = alignments[:i] + alignments[i+1:]

    for align1, align2, identity, gaps, score in alignments:
        length = len(align1)
        max_prefix_len = len(str(length)) + 1
        
        output_buffer('#===============================')
        output_buffer('# Parameters:')
        output_buffer('#   - alignment type: {}'.format(alignment['align_type']))
        output_buffer('#   - gap score: {}'.format(alignment['gap_score']))
        output_buffer('#   - mismatch score: {}'.format(alignment['mismatch_score']))
        output_buffer('#   - match score: {}'.format(alignment['match_score']))
        output_buffer('#')
        output_buffer('# Length: {}'.format(length))
        output_buffer('# Identity: {}/{} ({:.1f}%)'.format(identity, length, identity/length*100))
        output_buffer('# Gaps: {}/{} ({:.1f}%)'.format(gaps, length, gaps/length*100))
        output_buffer('# Score: {:.1f}'.format(score))
        output_buffer('#===============================')
        output_buffer('')
        for i in range(0, length, 60):
            prefix = str(i + 1)
            pl = len(prefix)
            l1 = prefix + ' ' * (max_prefix_len - pl)
            l2 = ' ' * max_prefix_len
            l3 = prefix + ' ' * (max_prefix_len - pl)
            end = min(i + 60, length)
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

seqs = []
for i in range(2):
    seq = input('Sequence {}: '.format(i + 1))
    if not seq or seq[0] == '>':
        seq = ''
    input_temp = 'ignored'
    while input_temp:
        input_temp = input()
        if input_temp and input_temp[0] != '>':
            seq += input_temp.upper().rstrip()
    seqs.append(seq)

gap_score = input('Gap score [-1]: ')
if gap_score:
    gap_score = float(gap_score)
else:
    gap_score = -1
mismatch_score = input('Mismatch score [0]: ')
if mismatch_score:
    mismatch_score = float(mismatch_score)
else:
    mismatch_score = 0
match_score = input('Match score [1]: ')
if match_score:
    match_score = float(match_score)
else:
    match_score = 1

alignment = build_matrix(seqs[0], seqs[1], gap_score, mismatch_score, match_score, align_type='local')
alignment = align_sequences(alignment, False)
print_alignment(alignment)
