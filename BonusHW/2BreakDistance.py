# The following script computes the 2-break distance between pair of genomes

import re
from collections import defaultdict

def get_component(node, g):
    q = [node]
    v = set()
    while q:
        node = q.pop(0)
        v.add(node)
        for n in g[node]:
            if n not in v:
                q += [n]
    return v

def read_genome_graph(s):
    g = defaultdict(list)
    for c in re.findall(r"\((.+?)\)", s):
        c = list(map(int, c.split()))
        for i in range(len(c) - 1):
            g[c[i]] += [-c[i + 1]]
            g[-c[i + 1]] += [c[i]]
        g[c[-1]] += [-c[0]]
        g[-c[0]] += [c[-1]]
    return g

# merge the graphs (assumes gp and gp have the same nodes)
def get_breakpoint_graph(p, q):
    g = {}
    for node in p.keys():
        g[node] = p[node] + q[node]
    return g

def get_break_distance(genomes):
    g = get_breakpoint_graph(*genomes)
    nodes = set(g.keys())
    blocks = len(nodes) // 2
    components = 0
    while len(nodes):
        res = get_component(next(iter(nodes)), g)
        nodes = nodes - res
        components += 1
    return blocks - components

def main(file):
    genomes = [read_genome_graph(s) for s in open(file).read().splitlines()]
    print(get_break_distance(genomes))

main("D:\dontm\Documents\School\B363\BonusHW\dataset.txt")