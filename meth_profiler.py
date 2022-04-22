import numpy as np
import random as random

def input_maker(methylations, sample_methylations, datasize, window_size, organism_name, half_w, threshold=0.5):
    mlevels = methylations['mlevel']
    mlevels = np.asarray(mlevels)
    chrs_counts = methylations['chr'].value_counts()
    last_chr_pos = {}
    chrnums = list(chrs_counts.index)
    sum = 0
    for i in range(len(chrnums)):
        if i in chrs_counts.keys():
            last_chr_pos[i] = sum+chrs_counts[i]-1
            sum += chrs_counts[i]
    # last_chr_pos ==> {0: 5524, 1: 1042784, 2: 1713034, 3: 2550983, 4: 3205486, 5: 4145381, 6: 4153872}
    # methylations.iloc[2550983] => chr 3.0 position    23459763.0
    # methylations.iloc[2550984] => chr 4.0 position    1007
    sub_methylations_me = sample_methylations[sample_methylations['mlevel'] >= threshold]
    sub_methylations_ume = sample_methylations[sample_methylations['mlevel'] < threshold]
    idxs_me = sub_methylations_me['idx']
    idxs_ume = sub_methylations_ume['idx']
    #if not check_indexes(methylations, sub_methylations_me) or not check_indexes(methylations, sub_methylations_ume):
    #    exit()
    X = np.zeros((datasize, window_size))
    Y = np.zeros(datasize)
    avlbls_me = np.asarray(idxs_me)
    avlbls_ume = np.asarray(idxs_ume)
    for lcp in list(last_chr_pos.values()):
        if lcp > 0 and lcp < len(mlevels) - window_size:
            avlbls_me = np.setdiff1d(avlbls_me, range(lcp-half_w, lcp+half_w))
            avlbls_ume = np.setdiff1d(avlbls_ume, range(lcp-half_w, lcp+half_w))
    #check_availables
    smple_me = random.sample(list(avlbls_me), int(datasize/2))
    smple_ume = random.sample(list(avlbls_ume), int(datasize/2))
    smple = np.concatenate([np.asarray(smple_me), np.asarray(smple_ume)])
    random.shuffle(smple)
    count_errored = 0
    #check_mlevels(methylations, mlevels)
    print('border conditions: ', np.count_nonzero(np.asarray(smple) < half_w))
    for index, p in enumerate(smple):
        try:
            X[index] = np.concatenate((mlevels[p-half_w: p], mlevels[p+1: p+half_w+1]), axis=0)
            Y[index] = 0 if mlevels[p] < 0.5 else 1
        except ValueError:
            #print(index, p)
            count_errored += 1
    X = X.reshape(list(X.shape) + [1])
    print(count_errored, ' profiles faced error')
    print('dataset balanceness', len(Y[Y > 0]), len(Y[Y==0]))
    return X, Y

def profiler(methylations, sample_methylations, datasize, organism_name, window_size=20):
    half_w = int(window_size/2)
    X, Y= input_maker(methylations, sample_methylations, int(datasize), window_size, organism_name, half_w)
    return X, Y
