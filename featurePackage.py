from functools import reduce

def Tfreq(seq):
    T_ratio = seq.count("T")/len(seq)
    return T_ratio

def codonFreq(seq):
    bases = ["A", "C", "G", "T"]
    codon64=[]
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon123 = base1 + base2 + base3
                codon64.append(codon123)
    phase_1 = range(0, len(seq), 3)
    codons = []
    freqList = []
    for i in phase_1:
        codon = seq[i:i+3]
        codons.append(codon)
    nums = len(codons)
    for i in codon64:
        freq = codons.count(i)/nums
        freqList.append(freq)
    return freqList

def singleZ(subSeq):
    fa = subSeq.count("A")/len(subSeq)
    fg = subSeq.count("G")/len(subSeq)
    fc = subSeq.count("C")/len(subSeq)
    ft = subSeq.count("T")/len(subSeq)
    x = (fa + fg) - (fc + ft)
    y = (fa + fc) - (fg + ft)
    z = (fa + ft) - (fg + fc)
    return x, y, z

def dinucleotideZ(subSeq1, subSeq2):
    bases = ["A", "C", "G", "T"]
    dibase_list = []
    dinucleotideZ_feature = []
    for i in range(0, len(subSeq1)):
        dibase_list.append(subSeq1[i]+subSeq2[i])
    dibase_list_len = len(dibase_list)
    for i in bases:
        XA = dibase_list.count(i + "A")/dibase_list_len
        XG = dibase_list.count(i + "G")/dibase_list_len
        XC = dibase_list.count(i + "C")/dibase_list_len
        XT = dibase_list.count(i + "T")/dibase_list_len
        XX = (XA + XG) - (XC + XT)
        YX = (XA + XC) - (XG + XT)
        ZX = (XA + XT) - (XG + XC)
        dinucleotideZ_feature.append(XX)
        dinucleotideZ_feature.append(YX)
        dinucleotideZ_feature.append(ZX)
    return dinucleotideZ_feature

def zcurve(sequence):
    phase_1 = range(0, len(sequence), 3)
    phase_2 = range(1, len(sequence), 3)
    phase_3 = range(2, len(sequence), 3)
    phase_1_seq = "".join([sequence[j] for j in phase_1])
    phase_2_seq = "".join([sequence[j] for j in phase_2])
    phase_3_seq = "".join([sequence[j] for j in phase_3])
    x1, y1, z1 = singleZ(phase_1_seq)
    x2, y2, z2 = singleZ(phase_2_seq)
    x3, y3, z3 = singleZ(phase_3_seq)
    dinucleotideZ_in12 = dinucleotideZ(phase_1_seq, phase_2_seq)
    dinucleotideZ_in23 = dinucleotideZ(phase_2_seq, phase_3_seq)
    dinucleotideZ_in13 = dinucleotideZ(phase_1_seq, phase_3_seq)
    zvariable = [x1, y1, z1, x2, y2, z2, x3, y3, z3]
    zvariable.extend(dinucleotideZ_in12)
    zvariable.extend(dinucleotideZ_in23)
    zvariable.extend(dinucleotideZ_in13)
    return zvariable

def countf(kstrings, k):
    kmers = reduce(lambda x, y: [i + j for i in x for j in y], [['A', 'T', 'C', 'G']] * k)
    strings_len = len(kstrings)
    res = []
    for i in kmers:
        i_frequency = round(kstrings.count(i)/strings_len, 4)
        res.append(i_frequency)
    return res

def kmer(sequence, k):
    sequence_len = len(sequence) - k + 1
    kstrings = []
    for i in range(0, sequence_len):
        kstrings.append(sequence[i:i+k])
    kmer_freq = countf(kstrings, k)
    return kmer_freq

