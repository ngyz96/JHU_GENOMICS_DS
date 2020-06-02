# What is the edit distance of the best match between pattern \verb|GCTGATCGATCGTACG|GCTGATCGATCGTACG and the excerpt of human chromosome 1? (Don't consider reverse complements.)
# What is the edit distance of the best match between pattern \verb|GATTTACCAGATTGAG|GATTTACCAGATTGAG and the excerpt of human chromosome 1? (Don't consider reverse complements.)

def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def editDistance_Fewest_Edits(P, T):
    # Create distance matrix
    D = []
    for i in range(len(P)+1):
        D.append([0]*(len(T)+1))
    # Initialize first column of matrix
    for i in range(len(P)+1):
        D[i][0] = i
    # Fill in the rest of the matrix
    for i in range(1, len(P)+1):
        for j in range(1, len(T)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if P[i-1] == T[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the min_value in the bottom row of the matrix
    # min_index = min(D[-1])
    return min(D[-1])

# def print_Matrix(D, P, T):
#     print('  ', end='')
#     for char in T:
#         print(char, end=' ')
#     print()
#     for i in range(len(D)):
#         for j in range(len(D[0])):
#             print(D[i][j], end=' ')
#         print('')

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

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlap_all_pairs(reads, k):
    #making kmer dic. key = kmer, value = set of reads with the particular kmer
    kmer_dic = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in kmer_dic.keys():
                kmer_dic[kmer] = set()
                kmer_dic[kmer].add(read)
            else:
                kmer_dic[kmer].add(read)
    #finding matches
    overlaps = []
    for read in reads:
        suffix = read[-k:]
        if suffix in kmer_dic.keys():
            set_reads = kmer_dic[suffix]
            for match in set_reads:
                val = overlap(read, match, k)
                if val > 0:
                    if read != match:
                        overlaps.append((read, match))
    return overlaps

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

# P = 'GCGTATGC'
# T = 'TATTGGCTATACGGTT'
# _, seq = readFasta('chr1.GRCh38.excerpt.fasta')
# T = seq[0]
# P = 'GCTGATCGATCGTACG'
# P = 'GATTTACCAGATTGAG'
# print(editDistance_Fewest_Edits(P,T))

# reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
# print(overlap_all_pairs(reads, 4))

sequences, _ = readFastq('ERR266411_1.for_asm.fastq')
overlaps = overlap_all_pairs(sequences, 30)
print(len(overlaps))
print(overlaps[:2])
suffix = set()
for overlap in overlaps:
        suffix.add(overlap[0])
print(len(suffix))