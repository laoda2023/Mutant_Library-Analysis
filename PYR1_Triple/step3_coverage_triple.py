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
            muts = tuple(muts.split(','))
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
    print('WT\t%5.3f %%' % (100 * wt / tot))
    print('Library encoded mutations:')
    for n in sorted(num):
        if n == 0:
            continue
        print('%i\t%5.3f %%' % (n, (100 * num[n] / tot)))
    print('Other mutations:\t%5.3f %%' % (100 * other_tot / tot))
    c_n_obs, c_n_tot, cov = coverage(target_muts, counts)
    print('Coverage of targeted triple mutations:\t%5.3f %% (%i / %i)' % (100*cov, c_n_obs, c_n_tot))
    print()
    #for truncation calculation
    with open('./pyr1_triple_count_q20.txt', 'r') as f:
        for _ in range(2):
            next(f)
        lines = f.readlines()
        col2 = [line.split()[1]for line in lines if '*' in line]
        col2_sum = sum(map(int,col2))
        print('Percentage with truncation:\t%5.3f' % (100*col2_sum/tot) , '% (other mutations includes truncation)') 
        
    print()
    print('------Summary of observed read counts------')
    print('total reads:%i' % (tot))# tot = lib_tot + other_tot
    print('total_reads_for_library:%i' % (lib_tot))
    print('total_reads_non_library:%i' % (other_tot))#it includes truncation
    print('WT_reads:%i' % (num[0]))
    print('Reads_with_single mutants:%i' % (num[1]))
    print('Reads_with_double mutants:%i' % (num[2]))
    print('Reads_with_triple mutants:%i' % (num[3]))
    print('Reads_with_quadruple mutants:%i' % (num[4]))
    print('Reads_with_quintuple mutants:%i' % (num[5]))
    print('Reads_with_sextuple_mutants:%i' % (num[6]))
    print('Reads_with_septuple_mutants:%i' % (num[7]))
    print('Reads_with_octuple_mutants:%i' % (num[8]))
    print('Reads_with_nonuple_mutants:%i' % (num[9]))
    print('Reads_with_truncation:%i' % col2_sum) 


def print_stats_combo(lib_muts, target_muts, counts):
    c_n_obs, c_n_tot, cov = coverage(target_muts, counts)
    return 100*cov, c_n_obs, c_n_tot


with open('./library_mut_triple.txt') as f:
    lib_muts = f.read().splitlines()


pyr1_triple = []
pyr1_triple_300 = []
pyr1_triple_030 = []
pyr1_triple_003 = []
pyr1_triple_111 = []
pyr1_triple_120 = []
pyr1_triple_102 = []
pyr1_triple_012 = []
pyr1_triple_021 = []
pyr1_triple_210 = []
pyr1_triple_201 = []

for mi in lib_muts:
    i = mut_pos(mi)
    for mj in lib_muts:
        j = mut_pos(mj)
        for mk in lib_muts:
        	k = mut_pos(mk)
        	if i == j or i == k or j == k:
        		continue
        	if i < j and j < k:
        		pyr1_triple.append((mi,mj,mk))
        		if i <= 94 and 108 <= j <= 141 and 159 <= k <= 167:
        			pyr1_triple_111.append((mi,mj,mk))
        		if 159 <= i <= 167 and 159 <= j <= 167 and 159 <= k <= 167:
        			pyr1_triple_003.append((mi,mj,mk)) 
        		if 108 <= i <= 141 and 108 <= j <= 141 and 108 <= k <= 141:
        			pyr1_triple_030.append((mi,mj,mk))
        		if i <= 94 and j <= 94 and k <= 94:
        			pyr1_triple_300.append((mi,mj,mk)) 
        		if i <= 94 and 108 <= j <= 141 and 108 <= k <= 141:
        			pyr1_triple_120.append((mi,mj,mk)) 
        		if i <= 94 and 159 <= j <= 167 and 159 <= k <= 167:
        			pyr1_triple_102.append((mi,mj,mk)) 
        		if i <= 94 and j <= 94 and 108 <= k <= 141:
        			pyr1_triple_210.append((mi,mj,mk)) 
        		if i <= 94 and j <= 94 and 159 <= k <= 167:
        			pyr1_triple_201.append((mi,mj,mk))
        		if 108 <= i <= 141 and 159 <= j <= 167 and 159 <= k <= 167:
        			pyr1_triple_012.append((mi,mj,mk)) 
        		if 108 <= i <= 141 and 108 <= j <= 141 and 159 <= k <= 167:
        			pyr1_triple_021.append((mi,mj,mk)) 
        		


pyr1_triple_kept, pyr1_triple_discarded, pyr1_triple_counts = \
    read_counts('counts_mutations_triple.txt')

print()
print('------Key statistical data------')
print_stats(lib_muts, pyr1_triple, pyr1_triple_counts)

print()
print('------Designed mutants for each combination------')
print('triple 1+1+1: ', len(pyr1_triple_111))
print('triple 0+0+3: ', len(pyr1_triple_003))
print('triple 0+3+0: ', len(pyr1_triple_030))
print('triple 3+0+0: ', len(pyr1_triple_300))
print('triple 1+2+0: ', len(pyr1_triple_120))
print('triple 1+0+2: ', len(pyr1_triple_102))
print('triple 2+1+0: ', len(pyr1_triple_210))
print('triple 2+0+1: ', len(pyr1_triple_201))
print('triple 0+1+2: ', len(pyr1_triple_012))
print('triple 0+2+1: ', len(pyr1_triple_021))
print('pyr1_triple:  ', len(pyr1_triple))

print()

print('------Coverage of each combination------')
print('triple 1+1+1: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_111, pyr1_triple_counts)))
print('triple 0+0+3: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_003, pyr1_triple_counts)))
print('triple 0+3+0: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_030, pyr1_triple_counts)))
print('triple 3+0+0: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_300, pyr1_triple_counts)))
print('triple 1+2+0: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_120, pyr1_triple_counts)))
print('triple 1+0+2: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_102, pyr1_triple_counts)))
print('triple 2+1+0: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_210, pyr1_triple_counts)))
print('triple 2+0+1: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_201, pyr1_triple_counts)))
print('triple 0+1+2: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_012, pyr1_triple_counts)))
print('triple 0+2+1: {:.3f}% ({}/{})'.format(*print_stats_combo(lib_muts, pyr1_triple_021, pyr1_triple_counts)))

