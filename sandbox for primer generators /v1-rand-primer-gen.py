import sys
import random
import multiprocessing

primer_library = []
MAX_HP = 2
PRIMER_BPS = 14
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6

complement_map = {'G': 'C',
                  'C': 'G',
                  'T': 'A',
                  'A': 'T'}


# here we try to take a different approach to the code by using set instead of list
# and each time we generate(random) a primer we check if it is in the set and "roll" off another one
# in addition we use parallel processing to speed up the process and generate more primers to reduce the time
# we use the same constants as michelle's code and 2-3 functions to use here
def complement_strand(strand):
    return ''.join(complement_map[base] for base in strand)


def cg_percent(strand):
    return (strand.count('C') + strand.count('G')) / len(strand)


def max_homopolymer(strand):
    max_hp = 1
    cur_hp = 1
    for i in range(1, len(strand)):
        if strand[i] == strand[i - 1]:
            cur_hp += 1
            max_hp = max(max_hp, cur_hp)
        else:
            cur_hp = 1
    return max_hp


def contains_complement(strand1, strand2, length):
    complement = complement_strand(strand1)
    for i in range(len(complement) - length):
        if complement[i: (i + length)] in strand2:
            return True
    return False


def hamming_distance(strand1, strand2):
    return sum(a != b for a, b in zip(strand1, strand2))


def next_primer(strand):  # different from michelle's code
    next_base = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
    return ''.join(next_base[base] for base in strand)


def generate_random_primers(num_samples):  # also different approach from michelle's code
    random.seed()
    primers = set()
    while len(primers) < num_samples:
        primer = ''.join(random.choice('ACGT') for _ in range(PRIMER_BPS))  # generate random primer of length 14
        if cg_percent(primer) >= 0.45 and cg_percent(primer) <= 0.55 and \
                max_homopolymer(primer) <= MAX_HP and \
                not contains_complement(primer, primer, MAX_SELF_COMP + 1) and \
                all(hamming_distance(primer, p) >= MIN_HAM for p in primer_library) and \
                all(not contains_complement(p, primer, MAX_INTER_COMP + 1) for p in primer_library):
            primers.add(primer)
    return primers


def add_primers_to_library(primer):
    i = 0;
    if all(hamming_distance(primer, p) >= MIN_HAM for p in primer_library) and \
            all(not contains_complement(p, primer, MAX_INTER_COMP + 1) for p in primer_library):
        i += 1
        primer_library.append(primer)
        if i % 50 == 0: sys.stdout.write("\n" + i + " primers sofar")


def generate_primers_parallel(num_samples, num_processes):
    pool = multiprocessing.Pool(processes=num_processes)
    results = pool.map(generate_random_primers, [num_samples // num_processes] * num_processes)
    for result in results:
        for primer in result:
            add_primers_to_library(primer)
    pool.close()
    pool.join()


if __name__ == "__main__":
    # maybe generate all 4^14 primers possible
    # minus the ones that are not valid for reasons like cg content, homopolymers, self complementing
    # and then when filter all the ones not possible by themselves from the set
    # and start  polling primers from the set in random order
    # and try getting max number of primers (hamming distance, inter strand complementing)

    generate_primers_parallel(1000000, 269)

    print("\nNumber of primers:", len(primer_library))
