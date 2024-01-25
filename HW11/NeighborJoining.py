from collections import defaultdict
import numpy as np

def parse_matrix(lines):
    """Parse an int matrix from a set of lines"""
    return [[int(x) for x in y.split()] for y in lines]

def get_edges(graph):
    e = []
    for i in sorted(graph):
        for j in graph[i]:
            e += [f"{i}->{j['n']}:{j['w']:.3f}"]
            e += [f"{j['n']}->{i}:{j['w']:.3f}"]
    return sorted(e)

def get_closest(D):
    D = np.copy(D)
    np.fill_diagonal(D, D.max() + 1)
    return divmod(D.argmin(), D.shape[1])

def joining_matrix(D, n):
    ND = np.copy(D)
    for i in range(len(D)):
        for j in range(len(D)):
            if i != j:
                ND[i, j] = (n - 2) * D[i, j] - sum(D[i, :]) - sum(D[j, :])
    return ND

def neighbor_joining(D, n, names=None):
    if not names:
        names = list(range(n))

    if n == 2:
        T = defaultdict(list)
        T[names[0]].append({"n": names[1], "w": D[0][1]})
        return T

    ND = joining_matrix(D, n)
    i, j = get_closest(ND)
    delta = (sum(D[i, :]) - sum(D[j, :])) / (n - 2)
    branch_i = (D[i, j] + delta) / 2
    branch_j = (D[i, j] - delta) / 2

    bi = names[i]
    bj = names[j]

    D = np.append(D, np.zeros((1, len(D))), axis=0)
    D = np.append(D, np.zeros((len(D), 1)), axis=1)
    names = names + [max(names) + 1]

    for k in range(n):
        D[k, n] = (D[k, i] + D[k, j] - D[i, j]) / 2
        D[n, k] = (D[k, i] + D[k, j] - D[i, j]) / 2
    for x in [j, i]:
        D = np.delete(D, x, 0)
        D = np.delete(D, x, 1)
        del names[x]

    T = neighbor_joining(D, n - 1, names)

    T[names[-1]].append({"n": bi, "w": branch_i})
    T[names[-1]].append({"n": bj, "w": branch_j})
    
    return T

def main(file):
    n, *D = open(file).read().splitlines()
    D = np.array(parse_matrix(D), float)
    g = neighbor_joining(D, int(n))
    for e in get_edges(g):
        print(e)

main("D:\dontm\Documents\School\B363\HW11\dataset.txt")