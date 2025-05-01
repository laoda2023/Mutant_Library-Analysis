def mut_pos(mut):
    return int(mut.split('-')[0][1:-1])

def read_counts(path):
    with open(path) as f:
        header, *lines = f.read().splitlines()
    _, kept, _, discarded = header.split()
    kept, discarded = int(kept), int(discarded)
    counts = {}
    for line in lines:
        muts, count = line.split('\t')
        count = int(count)
        if muts != '':
            muts = tuple(m.split('-')[0] for m in muts.split(','))#it deletes the codon
        if muts in counts:
            counts[muts] += count
        else:
            counts[muts] = count
    return kept, discarded, counts



def count_lib_only_muts(mutations, counts):
    total = sum(counts.values())
    num = {}
    for muts, k in counts.items():
        l = len(muts)
        if all(m in mutations for m in muts):
            num[l] = num.get(l, 0) + k
    lib_total = sum(num.values())
    other_total = total - lib_total
    return total, lib_total, other_total, num

def coverage(mutations, counts, min_count=1):
    n = len(mutations)
    k = len([m for m in mutations
             if (m in counts and counts[m] >= min_count)])
    return k, n, k/n

def print_stats(lib_muts, target_muts, counts):
    tot, lib_tot, other_tot, num = count_lib_only_muts(lib_muts, counts)
    wt = num[0]
    print('filtered reads\t%i' % sum(counts.values()))
    print('WT\t%5.2f' % (100 * wt / tot))
    print('Library encoded mutations:')
    for n in sorted(num):
        if n == 0:
            continue
        print('%i\t%5.2f' % (n, (100 * num[n] / tot)))
    print('Other mutations\t%5.2f' % (100 * other_tot / tot))
    c_n_obs, c_n_tot, cov = coverage(target_muts, counts)
    print('Coverage of targeted double mutations\t%5.2f (%i / %i)' % (100*cov, c_n_obs, c_n_tot))

with open('./library_mutations_DSM-Hao.txt') as f:#the style is K59S-TCT, containing codon.
    lib_muts = [line.split('-')[0] for line in f.read().splitlines()]


def print_stats_combo(lib_muts, target_muts, counts):
    c_n_obs, c_n_tot, cov = coverage(target_muts, counts)
    print('Coverage:\t%5.2f (%i / %i)' % (100*cov, c_n_obs, c_n_tot))

pyr1_double = []

pyr1_double_110 = []
pyr1_double_101 = []
pyr1_double_011 = []
pyr1_double_200 = []
pyr1_double_020 = []
pyr1_double_002 = []

           
           

for mi in lib_muts:
    i = mut_pos(mi)
    for mj in lib_muts:
        j = mut_pos(mj)
        if j == i:
            continue
        if i < j: 
            pyr1_double.append((mi,mj))
            if i <= 94 and 108 <= j <= 141:
                pyr1_double_110.append((mi,mj))
            if i <= 94 and 159 <= j <= 167:
                pyr1_double_101.append((mi,mj))
            if 108 <= i <= 141 and 159 <= j <= 167:
                pyr1_double_011.append((mi,mj))
            if i <= 94 and j <= 94:
                pyr1_double_200.append((mi,mj))
            if 108 <= i <= 141 and 108 <= j <= 141:
                pyr1_double_020.append((mi,mj))
            if 159 <= i <= 167 and 159 <= j <= 167:
                pyr1_double_002.append((mi,mj))
                

pyr1_double_dict = {
'x': pyr1_double,
'110': pyr1_double_110,
'101': pyr1_double_101,
'011': pyr1_double_011,
'200': pyr1_double_200,
'020': pyr1_double_020,
'002': pyr1_double_002,

}   
#insert a code for deleting predicted constitutive mutants,ok!-----------------
fp_predicted = 'predicted_1188fp.txt'

with open(fp_predicted) as f:
    for line_num, line in enumerate(f, 1):
        mutation = line.strip()
        if not mutation:
            continue
        target_mut = tuple(mutation.split('_'))

        if len(target_mut) != 2:
            print(f" the {line_num} row has style errors：'{mutation}' -> {target_mut}")
            continue  #if error occurred, interrupt
        
        #delete the matched constitutive mutants
        for key in pyr1_double_dict:
            if target_mut in pyr1_double_dict[key]:
                pyr1_double_dict[key].remove(target_mut)
#insert a code for deleting predicted constitutive mutants,ok!-----------------


print('pyr1_double: ', len(pyr1_double))
print('pyr1_double_110: ', len(pyr1_double_110))
print('pyr1_double_101: ', len(pyr1_double_101))
print('pyr1_double_011: ', len(pyr1_double_011))
print('pyr1_double_200: ', len(pyr1_double_200))
print('pyr1_double_020: ', len(pyr1_double_020))
print('pyr1_double_002: ', len(pyr1_double_002))

reads_kept, reads_discarded, reads_mut_counts = \
    read_counts('./counts_mut_DSM-Hao.txt') #the style is K59S-TCT, containing codon.


assert(reads_kept == sum(reads_mut_counts.values()))


print('pyr1_double:')
print_stats(lib_muts, pyr1_double, reads_mut_counts)

print('pyr1_double_110:')
print_stats_combo(lib_muts, pyr1_double_110, reads_mut_counts)
print('pyr1_double_101:')
print_stats_combo(lib_muts, pyr1_double_101, reads_mut_counts)
print('pyr1_double_011:')
print_stats_combo(lib_muts, pyr1_double_011, reads_mut_counts)
print('pyr1_double_200:')
print_stats_combo(lib_muts, pyr1_double_200, reads_mut_counts)
print('pyr1_double_020:')
print_stats_combo(lib_muts, pyr1_double_020, reads_mut_counts)
print('pyr1_double_002:')
print_stats_combo(lib_muts, pyr1_double_002, reads_mut_counts)



