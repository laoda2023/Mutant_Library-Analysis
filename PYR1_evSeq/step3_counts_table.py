import pandas as pd
from step1_merge import read_seqs, bytes_to_str
from tqdm import tqdm 

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

FIRST_CODON = 48
WT_ORF = 'ATCGTACGACGATTCGACAAACCACAAACATACAAACACTTCATCAAATCCTGCTCCGTCGAACAAAACTTCGAGATGCGCGTCGGATGCACGCGCGACGTGATCGTCATCAGTGGATTACCGGCGAACACATCAACGGAAAGACTCGATATACTCGACGACGAACGGAGAGTTACCGGATTCAGTATCATCGGAGGCGAACATAGGCTGACGAATTACAAATCCGTTACGACGGTGCATCGGTTCGAGAAAGAGAATCGGATCTGGACGGTGGTTTTGGAATCTTACGTCGTTGATATGCCGGAAGGTAACTCGGAGGATGATACTCGTATGTTTGCTGATACGGTTGTGAAGCTTAATTTGCAGAAACTCGCGACGGTTGCTGAAGCTATGGCTCGT'
WT_CODONS = cds_codons(WT_ORF)

primer_info = pd.read_csv("primer_map.csv")

result_df = pd.DataFrame(columns=["IndexPlate", "Well", "Reads_Counts", "Mutation_top1", "Count_top1", "Mutation_top2", "Count_top2", 
"Mutation_top3", "Count_top3",])

for index, primer_row in tqdm(primer_info.iterrows(), total = len(primer_info)): 
    fbc_seq = primer_row["FBC"]
    rbc_seq = primer_row["RBC"]
    fbc_name = primer_row["FBC-name"]
    rbc_name = primer_row["RBC-name"]
    plate_name = primer_row["IndexPlate"]
    Well_name = primer_row["Well"]
    matched_counts = 0
    mutation_counts = {}

    
    for seq_id, seq, qual_id, qual in read_seqs(open("./example_uniformed.fq")):
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
    
    top_mutations = [mutation[0] for mutation in sorted_mutation_counts[:3]]
    top_counts = [mutation[1] for mutation in sorted_mutation_counts[:3]]
 
    if len(sorted_mutation_counts) >= 3:
        new_row = pd.DataFrame([{
        "IndexPlate": plate_name,
        "Well": Well_name,
        "Reads_Counts": matched_counts,
        "Mutation_top1": top_mutations[0],
        "Count_top1": top_counts[0],
        "Mutation_top2": top_mutations[1],
        "Count_top2": top_counts[1],
        "Mutation_top3": top_mutations[2],
        "Count_top3": top_counts[2]
        }])
        result_df = pd.concat([result_df, new_row], ignore_index=True)


    elif len(sorted_mutation_counts) == 2:
        new_row = pd.DataFrame([{
        "IndexPlate": plate_name,
        "Well": Well_name,
        "Reads_Counts": matched_counts,
        "Mutation_top1": top_mutations[0],
        "Count_top1": top_counts[0],
        "Mutation_top2": top_mutations[1],
        "Count_top2": top_counts[1],
        "Mutation_top3": None,
        "Count_top3": None
        }])
        result_df = pd.concat([result_df, new_row], ignore_index=True)

    elif len(sorted_mutation_counts) == 1:
        new_row = pd.DataFrame([{
        "IndexPlate": plate_name,
        "Well": Well_name,
        "Reads_Counts": matched_counts,
        "Mutation_top1": top_mutations[0],
        "Count_top1": top_counts[0],
        "Mutation_top2": None,
        "Count_top2": None,
        "Mutation_top3": None,
        "Count_top3": None
        }])
        result_df = pd.concat([result_df, new_row], ignore_index=True)

    elif len(sorted_mutation_counts) == 0:
        new_row = pd.DataFrame([{
        "IndexPlate": plate_name,
        "Well": Well_name,
        "Reads_Counts": matched_counts,
        "Mutation_top1": None,
        "Count_top1": None,
        "Mutation_top2": None,
        "Count_top2": None,
        "Mutation_top3": None,
        "Count_top3": None
        }])
        result_df = pd.concat([result_df, new_row], ignore_index=True)
        
#save file
result_df.to_csv("mutations_results_top3.csv", index=False)
print("The results has been saved to mutations_results.csv")


