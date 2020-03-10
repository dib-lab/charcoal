#! /usr/bin/env python
import argparse
import numpy as np
from numpy import genfromtxt
import math
from scipy.cluster.hierarchy import dendrogram, linkage

from matplotlib import pyplot as plt
import pylab
from sourmash import fig


def load_and_normalize(filename):
    mat = genfromtxt(filename, delimiter=',')
    assert mat.shape[0] == 347            # number of metagenomes
    n_hashes = mat.shape[1]

    # go through and normalize all the sample-presence vectors for each hash;
    # track those with all 0s for later removal.
    to_delete = []
    for i in range(n_hashes):
        if sum(mat[:, i]):
            mat[:, i] /= math.sqrt(np.dot(mat[:, i], mat[:, i]))
        else:
            to_delete.append(i)

    # remove all columns with zeros
    print('removing {} null presence vectors'.format(len(to_delete)))
    for row_n in reversed(to_delete):
        mat = np.delete(mat, row_n, 1)

    assert mat.shape[1] == n_hashes - len(to_delete)
    n_hashes = mat.shape[1]

    # construct distance matrix of all sample-presence vectors x themselves.
    D = np.zeros((n_hashes, n_hashes))
    for i in range(n_hashes):
        for j in range(n_hashes):
            D[i][j] = np.dot(mat[:, i], mat[:, j])

    # done!
    return D
        
        
def cluster_and_plot(D, prefix):
    linked = linkage(D, 'single')

    plt.figure(figsize=(10, 7))
    Z = dendrogram(linked,
               orientation='top',
               distance_sort='descending',
               show_leaf_counts=True)

    idx1 = Z['leaves']
    print(len(idx1))
    D = D[idx1, :]
    D = D[:, idx1]

    im = plt.matshow(D, aspect='auto', origin='lower',
                           cmap=pylab.cm.YlGnBu)
    plt.savefig('{}.png'.format(prefix))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('matrix_csv')
    p.add_argument('output_fig')
    args = p.parse_args()

    mat = load_and_normalize(args.matrix_csv)

    n_hashes = mat.shape[0]
    labels = [""]*mat.shape[0]
    print('plotting {} hashes.'.format(n_hashes))

    x = fig.plot_composite_matrix(mat, labels,
                                  show_labels=False, show_indices=False,
                                  force=True)
    x.savefig(args.output_fig)


if __name__ == '__main__':
    main()
