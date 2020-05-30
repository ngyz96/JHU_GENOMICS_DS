from bm_preproc import BoyerMoore
from kmer_index import Index, SubseqIndex
def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    num_alignments = 0 
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        num_alignments += 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            num_character_comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons

def naive_with_counts(p, t):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        num_alignments += 1
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_alignments, num_character_comparisons

def readFasta(filename):
    meta_data = []
    sequences = []
    seq = ''
    with open(filename, 'r') as f:
        line = f.readline().rstrip()
        while True:
            if len(line) == 0:
                sequences.append(seq)
                break
            elif line[0] == '>':
                if seq:
                    sequences.append(seq)
                    seq = ''
                meta_data.append(line)
                line = f.readline().rstrip()
                continue
            else:
                seq = seq + line
                line = f.readline().rstrip()
                
    return meta_data, sequences

def naive_2mm(p, t):
    occurrences = []
    mm_counter = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mm_counter = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j] :  # compare characters
                mm_counter += 1
            if mm_counter > 2:
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def index_query(p, t, kmer):
    index = []
    q_len = 8
    q1,q2,q3 = p[:8], p[8:16], p[16:]
    hits1 = kmer.query(q1)
    hits2 = kmer.query(q2)
    hits3 = kmer.query(q3)

    for hit in hits1:
        if hit + len(p) > len(t):
            continue
        mismatch = 0 
        for i in range(q_len, len(p)):
            if p[i] != t[hit+i]:
                mismatch += 1
                if mismatch > 2:
                    break
        if mismatch <= 2:
            index.append((hit, hit + len(p)))

    for hit in hits2:
        if hit - q_len < 0 or hit + q_len > len(t):
            continue
        mismatch = 0
        for i in range(q_len):
            if p[i] != t[hit-q_len+i]:
                mismatch += 1
                if mismatch > 2:
                    break
        for j in range(2*q_len, len(p)):
            if p[j] != t[hit-q_len+j]:
                mismatch += 1
                if mismatch > 2:
                    break
        if mismatch <= 2:
            index.append((hit-q_len, hit + 2*q_len))
    for hit in hits3:
        if hit - 2*q_len < 0:
            continue
        mismatch = 0
        for i in range(2*q_len):
            if p[i] != t[hit-2*q_len+i]:
                mismatch += 1
            if mismatch > 2:
                break
        if mismatch <= 2:
            index.append((hit-2*q_len, hit + q_len))
    return index

# p = 'word'
# t = 'there would have been a time for such a word'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# p = 'needle'
# t = 'needle need noodle needle'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)

# p = 'word'
# t = 'there would have been a time for such a word'
# lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# p = 'needle'
# t = 'needle need noodle needle'
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)

# p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# _, seq = readFasta('chr1.GRCh38.excerpt.fasta')
# t = seq[0]
# nucleotides = 'ACGT'
# p_bm = BoyerMoore(p, nucleotides)
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)

# p = 'GGCGCGGTGGCTCACGCCTGTAAT'
# _, seq = readFasta('chr1.GRCh38.excerpt.fasta')
# t = seq[0]
# kmer = Index(t, 8)
# q1,q2,q3 = p[:8], p[8:16], p[16:]
# hits = kmer.query(q1)
# hits.extend(kmer.query(q2))
# hits.extend(kmer.query(q3))

# index = index_query(p, t, kmer)
# print(len(set(index)))

# q1,q2,q3 = p, p[1:], p[2:]
# sskmer = SubseqIndex(t,8,3)
# hits = sskmer.query(q1)
# hits.extend(sskmer.query(q2))
# hits.extend(sskmer.query(q3))
# print(len(hits))
t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print(occurrences) # [0, 14]
print(num_index_hits) #6

## TO DO:
def query_subseq(p, t, subseq_ind):
