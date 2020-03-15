#! /usr/bin/env python
"""
Plot togetherness.

Takes output of 'match_metagenomes.py', do presence-vector distance
calculations for hashes, and cluster/plot.
"""
import argparse

from matplotlib import pyplot as plt
import pylab
import scipy.cluster.hierarchy as sch
import collections

from utils import load_and_normalize


def plot_composite_matrix(D, labeltext, show_labels=True, show_indices=True,
                          vmax=1.0, vmin=0.0, force=False):
    """Build a composite plot showing dendrogram + distance matrix/heatmap.

    Returns a matplotlib figure."""
    if D.max() > 1.0 or D.min() < 0.0:
        print('This matrix doesn\'t look like a distance matrix - min value {}, max value {}'.format(D.min(), D.max()))
        if not force:
            raise ValueError("not a distance matrix")
        else:
            print('force is set; scaling to [0, 1]')
            D -= D.min()
            D /= D.max()

    if show_labels:
        show_indices = True

    fig = pylab.figure(figsize=(11, 8))
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])

    # plot dendrogram
    Y = sch.linkage(D, method='complete')

    dendrolabels = labeltext
    if not show_labels:
        dendrolabels = [str(i) for i in range(len(labeltext))]

    Z1 = sch.dendrogram(Y, orientation='left', labels=dendrolabels,
                        no_labels=not show_indices)
    ax1.set_xticks([])

    xstart = 0.45
    width = 0.45
    if not show_labels:
        xstart = 0.315
    scale_xstart = xstart + width + 0.01

    # plot matrix
    axmatrix = fig.add_axes([xstart, 0.1, width, 0.6])

    # (this reorders D by the clustering in Z1)
    idx1 = Z1['leaves']
    D = D[idx1, :]
    D = D[:, idx1]

    # show matrix
    im = axmatrix.matshow(D, aspect='auto', origin='lower',
                          cmap=pylab.cm.YlGnBu, vmin=vmin, vmax=vmax)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([scale_xstart, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)

    return fig


def augmented_dendrogram(*args, **kwargs):
    ddata = sch.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            plt.plot(x, y, 'ro')
            plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    return ddata

def annotated_dendro(mat):
    # this is what makes the distances
    Y = sch.linkage(mat, method='complete')

    fig = pylab.figure(figsize=(11, 8))

    Z = augmented_dendrogram(Y, orientation='top', no_labels=True)

    return fig, Z


def main():
    p = argparse.ArgumentParser()
    p.add_argument('matrix_csv')
    p.add_argument('output_fig')
    p.add_argument('--dendro', default=None)
    p.add_argument('--newick', default=None)
    args = p.parse_args()

    mat, n_orig_hashes = load_and_normalize(args.matrix_csv)

    n_hashes = mat.shape[0]
    labels = [""]*mat.shape[0] # could be loaded from .hashes file...
    print('plotting {} hashes.'.format(n_hashes))

    x = plot_composite_matrix(mat, labels,
                              show_labels=False, show_indices=False,
                              force=True)
    x.savefig(args.output_fig)

    if args.newick:
        Y = sch.linkage(mat, method='complete')

        ###

        rootnode, nodelist = sch.to_tree(Y, rd=True)

        def traverse(node, indent=' '):
            is_leaf = ' '
            if node.is_leaf():
                is_leaf = '*'
            print('XXX', indent, node.get_id(), is_leaf)
            if node.is_leaf():
                return '{}'.format(node.get_id())
            else:
                l = traverse(node.get_left(), indent=indent + ' ')
                r = traverse(node.get_right(), indent=indent + ' ')
                return "({},{}){}".format(l, r, node.get_id())

        with open(args.newick, 'wt') as fp:
            print(traverse(rootnode), file=fp)


    if args.dendro:
        y, Z = annotated_dendro(mat)
        y.savefig(args.dendro)

        CUT_POINT=2.0
        Y = sch.linkage(mat, method='complete')

        cluster_ids = sch.fcluster(Y, t=CUT_POINT, criterion='distance')
        Z = augmented_dendrogram(Y, orientation='top', no_labels=True)

        # now, get leaves and leaf labels
        idx1 = Z['leaves']
        new_labels = Z['ivl']

        # build clusters => sets of hashes.        
        clusters = collections.defaultdict(set)

        for i, k in enumerate(idx1):
            cluster_id = cluster_ids[k]
            clusters[cluster_id].add(new_labels[i])

        largest_n = 0
        for i in clusters:
            if len(clusters[i]) > largest_n:
                largest_n = len(clusters[i])
            #print('cluster {} is {} in size'.format(i, len(clusters[i])))

        print('Cut point: {}'.format(CUT_POINT))
        print('LARGEST cluster is {} of {} original hashes ({} nonempty)'.format(largest_n, n_orig_hashes, n_hashes))


if __name__ == '__main__':
    main()
