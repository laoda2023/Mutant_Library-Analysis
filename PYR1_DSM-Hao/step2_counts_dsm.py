import numpy as np
from step1_merge_dsm import read_seqs, bytes_to_str, reverse_complement, MIN_QUAL as FASTQ_MIN_QUAL

AA_CODONS = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'N': ('AAT', 'AAC'),
    'D': ('GAT', 'GAC'),
    'C': ('TGT', 'TGC'),
    'Q': ('CAA', 'CAG'),
    'E': ('GAA', 'GAG'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'H': ('CAT', 'CAC'),
    'I': ('ATT', 'ATC', 'ATA'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'K': ('AAA', 'AAG'),
    'M': ('ATG',),
    'F': ('TTT', 'TTC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    '*': ('TAA', 'TAG', 'TGA'),
}

TRANSLATE = {}

for aa, codons in AA_CODONS.items():
    for codon in codons:
        TRANSLATE[codon] = aa

def cds_codons(cds):
    assert len(cds) % 3 == 0
    return [cds[i:i+3] for i in range(0, len(cds), 3)]

#minimum quality of base
MIN_MIN_QUAL = 20

#length of wild type ORF
FRAME = 408

#the location of the first amino acid
FIRST_CODON = 47
WT_ORF = 'TCAATCGTACGACGATTCGACAAACCACAAACATACAAACACTTCATCAAATCCTGCTCCGTCGAACAAAACTTCGAGATGCGCGTCGGATGCACGCGCGACGTGATCGTCATCAGTGGATTACCGGCGAACACATCAACGGAAAGACTCGATATACTCGACGACGAACGGAGAGTTACCGGATTCAGTATCATCGGAGGCGAACATAGGCTGACGAATTACAAATCCGTTACGACGGTGCATCGGTTCGAGAAAGAGAATCGGATCTGGACGGTGGTTTTGGAATCTTACGTCGTTGATATGCCGGAAGGTAACTCGGAGGATGATACTCGTATGTTTGCTGATACGGTTGTGAAGCTTAATTTGCAGAAACTCGCGACGGTTGCTGAAGCTATGGCTCGTAACTCC'
positions_to_check = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160, 163, 164, 167]
positions_in_seq = [position - 47 for position in positions_to_check]
print('The position of amino acids in reads:', positions_in_seq)
final_positions = [pos * 3 + base for pos in positions_in_seq for base in [0, 1, 2]]
print('The position of bases in reads:', final_positions)

complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "M"}

WT_CODONS = cds_codons(WT_ORF)

def write_mutation_counts(path, outpath):
    discarded = 0
    kept_f = 0
    kept_r = 0
    mut_count = {}
    with open(path) as f:
        for seq_id, seq, qual_id, qual in read_seqs(f):
            seq = bytes_to_str(seq)
            qual_score = qual - FASTQ_MIN_QUAL
            
            if 'N' in seq:
                discarded += 1
                continue

            if discarded ==1 and 'TTACGAGC' in seq:
                #transform into the antiparallel sequences
                seq_r = "".join([complement[base] for base in seq[::-1]])
                #reverse the quality
                qual_r = qual[::-1]
				
                qual_list = []
                seq_list = []
                for pos in final_positions:
                    qual_list.append(qual_r[pos])
                    seq_list.append(seq_r[pos])
            #forward sequences
            if 'CAATCGTA' in seq and any(qual_score[pos] <=  MIN_MIN_QUAL for pos in final_positions):
                discarded += 1	   
                continue
                
            #forward sequences
            if 'CAATCGTA' in seq and all(qual_score[pos] >  MIN_MIN_QUAL for pos in final_positions):
			
                kept_f += 1   
                seq = seq[:FRAME]
                codons = cds_codons(seq)
                muts = tuple([(i, TRANSLATE[w], TRANSLATE[c], c)
					  for (i, (c, w)) in enumerate(zip(codons, WT_CODONS),
												   FIRST_CODON)
					  if c != w and TRANSLATE[c] != TRANSLATE[w]])
                mut_count[muts] = mut_count.get(muts, 0) + 1
                continue
			
			#reverse sequences
            if 'TTACGAGC' in seq:

                seq = "".join([complement[base] for base in seq[::-1]])
                qual = qual[::-1]
                qual = qual - FASTQ_MIN_QUAL
                
                if any(qual[pos] <=  MIN_MIN_QUAL for pos in final_positions):
			
                    discarded += 1
                    continue
			
                if all(qual[pos] >  MIN_MIN_QUAL for pos in final_positions):        
				
                    kept_r += 1   
			        
			        #remove the last base
                    seq = seq[:FRAME]
                    codons = cds_codons(seq)
                    muts = tuple([(i, TRANSLATE[w], TRANSLATE[c], c)
                          for (i, (c, w)) in enumerate(zip(codons, WT_CODONS),
                                                       FIRST_CODON)
                          if c != w and TRANSLATE[c] != TRANSLATE[w]])
                    mut_count[muts] = mut_count.get(muts, 0) + 1
#---------------------------------------------------
            	

            if (kept_f + kept_r + discarded) % 10000 == 0:
                print('All: ', kept_f + kept_r + discarded)
                print("kept_r: ", kept_r)
                print("kept_f: ", kept_f)
                print("discarded: ", discarded)      	

    with open(outpath, 'w') as f:
        print('kept_f_r %i discarded %i' % (kept_f + kept_r, discarded), file=f)
 
        for (muts, count) in sorted(mut_count.items(),
                                    key=lambda x: x[1],
                                    reverse=True):
            #If you want the generated file to exclude codons, use the following line of code.
            #s = "_".join([f"{wt}{pos}{mut}" for pos, wt, mut, codon in muts])
            
            #If you want the generated file to include codons, use the following line of code.
            s = ','.join(['%s%i%s-%s' % (w, i, m, c) for (i, w, m, c) in muts])
            print('\t'.join([s, str(count)]), file=f)


fastq_file_merged = './DSM_Hao_merged.fq'
write_mutation_counts(fastq_file_merged, 'counts_mutations_DSM-Hao.txt')