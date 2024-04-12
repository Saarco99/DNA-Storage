import itertools
import random
import multiprocessing

primer_library = set()
MAX_HP = 2
PRIMER_BPS = 14
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6

complement_map = {'G': 'C',
                  'C': 'G',
                  'T': 'A',
                  'A': 'T'}

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


def is_valid(primer):
    return cg_percent(primer) >= 0.45 and cg_percent(primer) <= 0.55 and \
           max_homopolymer(primer) <= MAX_HP and \
           not contains_complement(primer, primer, MAX_SELF_COMP + 1)


def poll_primers(valid_primers):#todo  consider adding parallelism here !!1
    primer_list = list(valid_primers)
    random.shuffle(primer_list)
    for primer in primer_list:
        if all(hamming_distance(primer, p) >= MIN_HAM for p in primer_library) and \
           all(not contains_complement(p, primer, MAX_INTER_COMP + 1) for p in primer_library):
            primer_library.add(primer)

# def worker(num_primers):
#     nucleotides = ['A', 'C', 'G', 'T']
#     all_primers = [''.join(p) for p in itertools.product(nucleotides, repeat=PRIMER_BPS)]
#     valid_primers = set(primer for primer in all_primers if is_valid(primer))
#     return valid_primers
def worker(num_primers):#with progress tracking
    nucleotides = ['A', 'C', 'G', 'T']
    all_primers = [''.join(p) for p in itertools.product(nucleotides, repeat=PRIMER_BPS)]
    valid_primers = set()
    quarter = len(all_primers) // 4
    for i, primer in enumerate(all_primers):
        if is_valid(primer):
            valid_primers.add(primer)
        if (i + 1) % quarter == 0:
            print(f"Worker {multiprocessing.current_process().name} has processed {(i + 1) / len(all_primers) * 100}% of its primers")
    return valid_primers

def generate_all_possible_primers_parallel(num_processes):
    pool = multiprocessing.Pool(processes=num_processes)
    num_primers_per_process = 4**PRIMER_BPS // num_processes
    results = pool.map(worker, [num_primers_per_process] * num_processes)
    pool.close()
    pool.join()
    return set().union(*results)

if __name__ == "__main__":
    num_processes = multiprocessing.cpu_count()
    valid_primers = generate_all_possible_primers_parallel(num_processes-4)
    poll_primers(valid_primers)
    print("\nNumber of primers:", len(primer_library))