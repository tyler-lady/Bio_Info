from random import randint, choices

import numpy

def random_select_motifs(Dna, k, t):
  Motifs = []
  for seq in Dna:
    index = randint(0, len(seq) - k)
    Motifs.append(seq[index:index+k])

  return Motifs

def Profile(Motifs, k):
  profile = []
  for i in range(k):
    for j in range(len(Motifs)):
      if j == 0:
        profile.append({ 'A': 1, 'T': 1, 'C': 1, 'G': 1 })
      profile[i][Motifs[j][i]] += 1
  return profile

def motif(profile, Seq, k):
    nuc_loc = {nucleotide: index for index, nucleotide in enumerate('ACGT')}
    probs = []
    for i in range(len(Seq) - k):
        current_prob = 1.
        for j, nucleotide in enumerate(Seq[i:i + k]):
            current_prob *= profile[j][nuc_loc[nucleotide]]
        probs.append(current_prob)

    i = numpy.random.choice(len(probs), p = numpy.array(probs) / numpy.sum(probs))
    return Seq[i:i + k]

def score(Motifs, k, t):
  profile = Profile(Motifs, k)
  score = 0
  for a in range(len(profile)):
    score += (4 + t - profile[a][max(profile[a], key=profile[a].get)])
    # print(profile[a][max(profile[a])], 4 + t, score)
  return score

def profile_with_pseudocounts(motifs):
    prof = []
    for i in range(len(motifs[0])):
        col = ''.join([motifs[j][i] for j in range(len(motifs))])
        prof.append([float(col.count(nuc)+1)/float(len(col)+4) for nuc in 'ACGT'])
    return prof

def gibbs_sampler(Dna, k, t, N):
  motifs = random_select_motifs(Dna, k, t)
  BestMotifs = list(motifs)
  for j in range(N):
    i = randint(0, t - 1)
    profile = profile_with_pseudocounts([x for amotif, x in enumerate(motifs) if amotif != i])#Profile(Motifs, k)
    motifs[i] = motif(profile, Dna[i], k)

    if score(motifs, k, t) < score(BestMotifs, k, t):
      BestMotifs = list(motifs)
  return BestMotifs

with open('D:\Downloads\\rosalind_ba2g.txt', 'r') as f:
  k, t, N = [int(x) for x in f.readline().strip().split(' ')]
  Dna = [x.strip() for x in f.readlines()]

  BestMotifs = gibbs_sampler(Dna, k, t, N)

  for i in range(0, 20):
    Motifs = gibbs_sampler(Dna, k, t, N)
    if score(Motifs, k, t) < score(BestMotifs, k, t):
      BestMotifs = Motifs
      print(BestMotifs, score(BestMotifs,k, t))
with open('rosalind_GS_output.txt', 'w') as out:
  for motif in BestMotifs:
    out.write(motif + '\n')