# The following script implements hierarchical clustering

import numpy as np
from collections import defaultdict

# find first min off diagonal index in array
def get_closest(D):
    D = np.copy(D)
    np.fill_diagonal(D, D.max() + 1)
    return divmod(D.argmin(), D.shape[1])

# replace the ith row/col with the average of the ith and jth and remove the jth
def average_replacement(D, i, j, di, dj):
    D = np.copy(D)
    av = (D[i, :] * di + D[j, :] * dj) / (di + dj)
    D[i, :] = av
    D[:, i] = av
    D = np.delete(D, j, 0)
    D = np.delete(D, j, 1)
    np.fill_diagonal(D, 0)
    return D

def get_children(T, node):
    q = [node]
    x = []
    while len(q):
        n = q.pop(0)
        if n in T:
            q += T[n]
        else:
            x += [n]
    return x


def hierarchical_clustering_generator(D, n):
    clusters = list(range(1, n + 1))
    T = {}
    size = defaultdict(lambda: 1)  # the number of children of a node
    node = n
    while len(clusters) > 1:
        node += 1
        i, j = get_closest(D)
        a, b = clusters[i], clusters[j]
        T[node] = [a, b]
        size[node] = size[a] + size[b]
        D = average_replacement(D, *get_closest(D), size[a], size[b])
        clusters[i] = node
        del clusters[j]
        yield get_children(T, a) + get_children(T, b)


def main(file):
    n, *m = open(file).read().splitlines()
    m = np.array([list(map(float, x.split())) for x in m])
    for step in hierarchical_clustering_generator(m, int(n)):
        print(*step)

main("D:\dontm\Documents\School\B363\BonusHW\dataset.txt")