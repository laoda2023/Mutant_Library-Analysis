import pandas as pd
import multiprocessing
from step1_merge import read_seqs, bytes_to_str
from tqdm import tqdm  # for progress of program

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

FIRST_CODON = 48 # I48
WT_ORF = 'ATCGTACGACGATTCGACAAACCACAAACATACAAACACTTCATCAAATCCTGCTCCGTCGAACAAAACTTCGAGATGCGCGTCGGATGCACGCGCGACGTGATCGTCATCAGTGGATTACCGGCGAACACATCAACGGAAAGACTCGATATACTCGACGACGAACGGAGAGTTACCGGATTCAGTATCATCGGAGGCGAACATAGGCTGACGAATTACAAATCCGTTACGACGGTGCATCGGTTCGAGAAAGAGAATCGGATCTGGACGGTGGTTTTGGAATCTTACGTCGTTGATATGCCGGAAGGTAACTCGGAGGATGATACTCGTATGTTTGCTGATACGGTTGTGAAGCTTAATTTGCAGAAACTCGCGACGGTTGCTGAAGCTATGGCTCGT'
WT_CODONS = cds_codons(WT_ORF)

#read the table with barcoded primers
primer_info = pd.read_csv("./primer_map.csv")

result_df = pd.DataFrame(columns=["IndexPlate", "Well", "Reads_Counts"])
all_column_names = []
all_mutations = []

for index, primer_row in tqdm(primer_info.iterrows(), total = len(primer_info)):
    fbc_seq = primer_row["FBC"]
    rbc_seq = primer_row["RBC"]
    fbc_name = primer_row["FBC-name"]
    rbc_name = primer_row["RBC-name"]
    plate_name = primer_row["IndexPlate"]
    Well_name = primer_row["Well"]
    matched_counts = 0
    mutation_counts = {}
    sorted_mutation_counts = []

    for seq_id, seq, qual_id, qual in read_seqs(open("example1_uniformed.fq")):
        seq = bytes_to_str(seq)
        if seq.startswith(fbc_seq) and seq.endswith(rbc_seq):
            matched_counts += 1
            seq = seq[9:-9]
            codons = cds_codons(seq)
            muts = tuple([(i, TRANSLATE[w], TRANSLATE[c], c)
                for (i, (c, w)) in enumerate(zip(codons, WT_CODONS),
                                                       FIRST_CODON)
                if c != w and TRANSLATE[c] != TRANSLATE[w]])

            mutation_name = "_".join([f"{wt}{pos}{mut}" for pos, wt, mut, codon in muts])
            mutation_counts[mutation_name] = mutation_counts.get(mutation_name, 0) + 1
            sorted_mutation_counts = sorted(mutation_counts.items(), key = lambda x: x[1], reverse = True)
    
    top_mutations = [mutation[0] for mutation in sorted_mutation_counts]
    top_counts = [mutation[1] for mutation in sorted_mutation_counts]
    top_mutations = ["WT" if mutation.strip() == '' else mutation for mutation in top_mutations]#if no mutation, put "WT" in the table.
    
    #require the counts of sorted mutations
    max_mutation_column = len(sorted_mutation_counts)
    mutation_and_count_names = [f"Mutation_top{i+1}" for i in range(max_mutation_column)] + [f"Count_top{i+1}" for i in range(max_mutation_column)]
    all_column_names.extend(mutation_and_count_names)
    
    result_dict = {
        "IndexPlate": plate_name,
        "Well": Well_name,
        "Reads_Counts": matched_counts,
    }

    for i, mutation in enumerate(top_mutations):
        result_dict[f"Mutation_top{i+1}"] = mutation
        result_dict[f"Count_top{i+1}"] = top_counts[i]
    
    #result_df = result_df.append(result_dict, ignore_index=True)
    result_df = pd.concat([result_df, pd.DataFrame([result_dict])], ignore_index=True)
    max_mutation_column = len(sorted_mutation_counts)
    all_mutations.append(max_mutation_column)
    
max_column = max(all_mutations) #find the max column

#generate all the columns
mutation_and_count_names = [f"Mutation_top{i+1}" for i in range(max_column)] + [f"Count_top{i+1}" for i in range(max_column)]
mutation_and_count_names = [name for pair in zip(mutation_and_count_names[:max_column], mutation_and_count_names[max_column:]) for name in pair]

#add other columns into generated columns
result_df = result_df.reindex(columns=["IndexPlate", "Well", "Reads_Counts"] + mutation_and_count_names)

#save results into files
result_df.to_csv("mutations_results_top_all.csv", index=False)
print("All the mutations have been saved into mutations_results.csv")