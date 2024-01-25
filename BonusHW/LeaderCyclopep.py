# The following script implements leaderboard cyclopeptide sequencing

import yaml

def get_mass():
    path = "mass.yaml"
    with open(path) as stream:
        return yaml.safe_load(stream)

def get_cyclo_spectrum(pep):
    c_spec = [0, sum(pep)]
    for i in range(1, len(pep)):
        for j in range(len(pep)):
            c_spec += [sum((pep[j:] + pep[:j])[0:i])]
    return sorted(c_spec)

def generate_substrings(text, size):
    for i in range(len(text) - size + 1):
        yield text[i : i + size]

def get_score(theoretical, expected):
    spec, score = theoretical[:], 0
    if spec:
        for m in expected:
            if m in spec:
                score += 1
                spec.remove(m)
    return score

def expand(peps, masses):
    return [p + [x] for x in masses for p in list(peps)]

def get_linear_spectrum(pep):
    spec = [0]
    for i in range(1, len(pep) + 1):
        for x in generate_substrings(pep, i):
            spec.append(sum(x))
    return spec

def get_linear_score(pep, spec):
    return get_score(get_linear_spectrum(pep), spec)

def get_cyclo_score(pep, spec):
    return get_score(get_cyclo_spectrum(pep), spec)

def cut(peps, spec, n):
    if len(peps) < n:
        return peps
    score = [get_linear_score(p, spec) for p in peps]
    lim = sorted(score, reverse=True)[n - 1]
    return [p for p, score in zip(peps, score) if score >= lim]

def get_leaderboard_cyclopeptide_sequencing(spec, n, masses):
    lb = [[]]
    lead = []
    while len(lb):
        lb = expand(lb, masses)
        for pep in lb.copy():
            if sum(pep) == spec[-1]:
                if get_cyclo_score(pep, spec) > get_cyclo_score(lead, spec):
                    lead = pep
            elif sum(pep) > spec[-1]:
                lb.remove(pep)
        lb = cut(lb, spec, n)
    return lead

def main(file):
    n, m = open(file).read().splitlines()
    spectrum = list(map(int, m.split()))
    masses = set(get_mass().values())
    p = get_leaderboard_cyclopeptide_sequencing(spectrum, int(n), masses)
    print(*p, sep="-")

main("D:\dontm\Documents\School\B363\BonusHW\dataset.txt")