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

def obtainIdentifiers(meta_data):
    identifiers = []
    for data in meta_data:
        identifier = ''
        data = data[1:]
        while data[0] != ' ':
            identifier += data[0]
            data = data[1:]
        identifiers.append(identifier)
    return identifiers

def findORF(sequences, frame):
    ORF = []
    length = []
    start_codon = 'ATG'
    stop_codon = ['TAA', 'TGA', 'TAG']
    for seq in sequences:
        for i in range(frame-1, len(seq) - 3 + 1, 3):
            if seq[i:i+3] == start_codon:
                for j in range(i+3, len(seq), 3):
                    if seq[j:j+3] in stop_codon:
                        ORF.append((i, j+2))
                        length.append(j+3-i)
                        break
    return ORF, length
def ComputingFrequencies(sequences, k):
    FrequencyArray = {}

    for seq in sequences:
        for i in range(len(seq) - k + 1):
            pattern = seq[i:i + k]
            if pattern not in FrequencyArray.keys():
                FrequencyArray[pattern] = 1
            else:
                FrequencyArray[pattern] += 1

    return FrequencyArray

meta_data, sequences = readFasta('dna2.fasta')
# Q1
# print(len(sequences))

# Q2
# length = []
# for seq in sequences:
#     length.append(len(seq))
# print(max(length))

# Q3
# print(min(length))

# Q4
# ORF, lengths = findORF(sequences, 2)
# print(max(lengths))

# Q5 
# ORF, lengths = findORF(sequences, 3)
# max_length = max(lengths)
# index = []
# for i, length in enumerate(lengths):
#     if length == max_length:
#         index.append(i)
# for j in index:
#     print(ORF[j][0] + 1)

# Q6
# _, ORF1_len = findORF(sequences, 1)
# _, ORF2_len = findORF(sequences, 2)
# _, ORF3_len = findORF(sequences, 3)

# ORF1_len.extend(ORF2_len)
# ORF1_len.extend(ORF3_len)
# print(max(ORF1_len))

# Q7
# identifiers = obtainIdentifiers(meta_data)
# index = None
# for i, identifier in enumerate(identifiers):
#     if identifier == 'gi|142022655|gb|EQ086233.1|16':
#         index = i

# _, ORF1_len = findORF([sequences[index]], 1)
# _, ORF2_len = findORF([sequences[index]], 2)
# _, ORF3_len = findORF([sequences[index]], 3)

# ORF1_len.extend(ORF2_len)
# ORF1_len.extend(ORF3_len)
# print(max(ORF1_len))

# Q8
# freq_array = ComputingFrequencies(sequences, 6)
# max_freq = max(freq_array.values())
# print(max_freq)

# Q9
# freq_array = ComputingFrequencies(sequences, 12)
# max_freq = max(freq_array.values())
# counter = 0
# for value in freq_array.values():
#     if value == max_freq:
#         counter += 1
# print(counter)

# Q10
# freq_array = ComputingFrequencies(sequences, 7)
# options = ['CATCGCC', 'CGCGCCG', 'GCGGCCG', 'TGCGCGC']
# for option in options:
#     if option in freq_array.keys():
#         print(option, freq_array[option])