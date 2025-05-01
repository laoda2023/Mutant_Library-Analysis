import numpy as np
from step1_merge import read_seqs, bytes_to_str, reverse_complement, MIN_QUAL as FASTQ_MIN_QUAL

MIN_MIN_QUAL = 20
complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

path = './example1_merged.fq'

discarded = 0
kept_f = 0
kept_r = 0
mut_count = {}
with open(path) as f:
    for seq_id, seq, qual_id, qual in read_seqs(f):
        seq = bytes_to_str(seq)
        qual_score = qual - FASTQ_MIN_QUAL
        if qual_score.min() < MIN_MIN_QUAL:
            discarded += 1
            continue 
        qual = bytes_to_str(qual)
        if 'CAATCGTA' in seq:
            kept_f += 1
            with open('example1_merged_f.fq', 'a') as f:
                f.write(seq_id + '\n')
                f.write(seq + '\n')
                f.write(qual_id + '\n')
                f.write(qual + '\n')

        if 'TTACGAGC' in seq:
            kept_r += 1
            seq = "".join([complement[base] for base in seq[::-1]])
            qual = qual[::-1]
        	
            with open('example1_merged_r.fq', 'a') as f:
                f.write(seq_id + '\n')
                f.write(seq + '\n')
                f.write(qual_id + '\n')
                f.write(qual + '\n')                  
                 	
        if (kept_f + kept_r + discarded) % 10000 == 0:
            print(kept_f + kept_r + discarded)
            print(kept_f)
            print(kept_r)
            print(discarded)