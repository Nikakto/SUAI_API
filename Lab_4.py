from collections import defaultdict, Counter, OrderedDict
import matplotlib.pyplot as plt

KMER_LENGTH = 25
TRASH = 5


if __name__ == '__main__':

    with open('hw3_dataset.fasta', 'r') as file:
        reads = file.read()

    reads = reads.split('\n')
    reads = reads[1::2]

    # ==================================================================================================================
    # DEBUG
    # ==================================================================================================================
    # KMER_LENGTH = 3
    # TRASH = 3
    # read = 'ATGGCTAGATCCCTACG'
    # reads = [ read[i:i+5] for i in range(len(read) - 5 + 1)]

    kmers = [ s[i:i+KMER_LENGTH] for s in reads for i in range(len(s) - KMER_LENGTH + 1) ]

    _counter = Counter(kmers)
    _unique = _counter.keys()
    _frequency = _counter.values()

    kmers_stat = dict(zip(_unique, _frequency))

    _counter = Counter(_frequency)
    kmers_frequency_data = dict(zip(_counter.keys(), _counter.values()))
    kmers_frequency = defaultdict(int, kmers_frequency_data)

    plt.figure()
    _max_x = max(kmers_frequency.keys())
    _x = list(range(_max_x+1))
    _y = [kmers_frequency[i] * i for i in range(_max_x+1)]
    plt.plot(_x, _y)
    plt.title(r'Assembly Dublication (k-mer length = %s)' % KMER_LENGTH)
    plt.xlabel('kmer multiplicity')
    plt.ylabel('Distinct Kmer Count')
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

    kmers_frequency = OrderedDict(sorted(kmers_frequency.items()))
    coverage, coverage_count = max(list(kmers_frequency.items())[TRASH:], key=lambda x: x[1])
    print('K-mer coverage:', coverage)

    genom_length = sum(kmers_frequency[i] * i for i in range(TRASH, len(kmers_frequency))) / coverage
    print('Approximate length of genome:', round(genom_length))