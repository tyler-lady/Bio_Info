# The following script reconstructs a string from its paired composition

from collections import defaultdict
from copy import deepcopy

def get_genome_path(dna):
    return dna[0] + "".join(x[-1] for x in dna[1:])

def find_cycle(graph, key):
    c = []
    c += [key]
    while len(graph[key]):
        key = graph[key].pop(0)
        c += [key]
    for k in list(graph):
        if len(graph[k]) == 0:
            del graph[k]
    return c

# Find an Eulerian cycle (if it exists).
def find_eulerian_cycle(graph):
    g = deepcopy(graph)
    c = find_cycle(g, list(g.keys())[0])
    while len(g):
        key = [x for x in c if x in g][0]
        i = c.index(key)
        c = c[:i] + find_cycle(g, key) + c[(i + 1) :]
    return c

def count_connections(graph):
    counts = defaultdict(lambda: {"in": 0, "out": 0})
    v = sum(graph.values(), [])
    for n in sorted(set(v)):
        counts[n]["in"] = v.count(n)
    for n in graph:
        counts[n]["out"] = len(graph[n])
    return counts

def get_eulerian_path(graph):
    graph = deepcopy(graph)
    conns = count_connections(graph)
    i = [k for k in conns if conns[k]["in"] > conns[k]["out"]][0]
    j = [k for k in conns if conns[k]["in"] < conns[k]["out"]][0]
    graph[i] = [j]
    cycle = find_eulerian_cycle(graph)[:-1]
    index = cycle.index(j)
    return cycle[index:] + cycle[:index]

def get_dbrujin_paired(pairs):
    g = defaultdict(list)
    for x in pairs:
        p = tuple([x[0][:-1], x[1][:-1]])
        s = tuple([x[0][1:], x[1][1:]])
        g[p].append(s)
    return g

def get_string_from_paired_composition(pairs, k, d):
    path = get_eulerian_path(get_dbrujin_paired(pairs))
    a = get_genome_path([x[0] for x in path])
    b = get_genome_path([x[1] for x in path])
    return a + b[-(k + d) :]

def main(file):
    ints, *pairs = open(file).read().splitlines()
    k, d = map(int, ints.split())
    pairs = [x.split("|") for x in pairs]
    print(get_string_from_paired_composition(pairs, k, d))

main("D:\dontm\Documents\School\B363\BonusHW\dataset.txt")