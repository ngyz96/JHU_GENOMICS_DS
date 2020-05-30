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

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive_with_rc(p, t):
    rc = reverseComplement(p)
    p_len = len(p)
    occurrences = []
    for i in range(len(t) - p_len + 1):  # loop over alignments
        if t[i:i+p_len] == p:
            occurrences.append(i)
        elif t[i:i+p_len] == rc:
            occurrences.append(i)

    return occurrences

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

def phred33ToQ(qual):
    return ord(qual) - 33

genome = readGenome('lambda_virus.fa')
occurrences = naive_2mm('AGGAGGTT', genome)
print(min(occurrences))

###### naive_with_rc test cases
# p = 'CCC'
# ten_as = 'AAAAAAAAAA'
# t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
# occurrences = naive_with_rc(p, t)
# print(occurrences)
# p = 'CGCG'
# t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
# occurrences = naive_with_rc(p, t)
# print(occurrences)
# phix_genome = readGenome('phix.fa')
# occurrences = naive_with_rc('ATTA', phix_genome)
# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))

###### naive_2mm test cases
# p = 'CTGT'
# ten_as = 'AAAAAAAAAA'
# t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
# occurrences = naive_2mm(p, t)
# print(occurrences)
# phix_genome = readGenome('phix.fa')
# occurrences = naive_2mm('GATTACA', phix_genome)
# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))



# _ , qualities = readFastq("ERR037900_1.first1000.fastq")
# pos = {}
# for x in qualities:
#     for i in range(len(x)):
#         quality = phred33ToQ(x[i])
#         if quality not in pos.keys():
#             pos[quality] = [i]
#         else:
#             pos[quality].append(i)
# for key, value in pos.items():
#     print(key, len(value))
# print(pos.keys())
# num = {}
# for x in pos[9]:
#     if x not in num.keys():
#         num[x] = 1
#     else:
#         num[x] += 1
# for key, value in num.items():
#     print(key, value)
# pos = [0]*100
# for q in qualities:
#     for i in range(len(q)):
#         pos[i] += phred33ToQ(q[i])
# import matplotlib.pyplot as plt
# plt.plot(range(len(pos)), pos)
# plt.show()
