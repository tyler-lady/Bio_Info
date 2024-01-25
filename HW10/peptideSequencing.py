from math import sqrt

def get_protein_weight_dict():
	table ='''A   71.03711
	C   103.00919
	D   115.02694
	E   129.04259
	F   147.06841
	G   57.02146
	H   137.05891
	I   113.08406
	K   128.09496
	L   113.08406
	M   131.04049
	N   114.04293
	P   97.05276
	Q   128.05858
	R   156.10111
	S   87.03203
	T   101.04768
	V   99.06841
	W   186.07931
	Y   163.06333''' 

	weight_dict = dict()

	for p in table.split('\n'):
		weight_dict[p.strip('\t').split()[0]] = float(p.strip('\t').split()[1])

	return weight_dict

# returns list with all the words possible from add list using add chars
def append_c(add_list, add_chars):
	newlist = []
	for item in add_list:
		newlist += [item+ch for ch in set(add_chars)]
	return newlist

# returns spectrum of peptide
def get_cyclic_spectrum(pep):
	# get weight dict
	w = get_protein_weight_dict()
	# init as mass 0 and mass of peptide
	spec = [0, sum([int(w[p]) for p in pep])]
	# find masses of adjacent subpeptides
	spec += [sum([int(weights[p]) for p in pep[j:j+i]]) for i in range(1,len(pep)) for j in range(len(pep)-i+1)]
	# sort list ascending order and change to strings
	spec = map(str,sorted(spec))

	return spec

with open('dataset.txt') as input:
	cyclospec = input.read().strip().split()

# get protein weight dict
weights = get_protein_weight_dict()

# n = length of peptide L = length of cyclospectrum (L = n(n-1) + 2)
# quadratic formula for n:  n = (sqrt(4L-7) + 1)/2
n = int((sqrt(4*len(cyclospec)-7)+1)/2)

# find first n protein in peptide  
prot, index = [], 1
while len(prot) != n:
	if int(cyclospec[index]) in map(int,weights.values()):
		prot.append(cyclospec[index])
	index += 1

# get names of proteins w/ given weight
names = []
for w in prot:
	names.append([items[0] for items in weights.items() if int(items[1])==int(w)][0])

# get possible seq
seq = append_c(names,names)
for rep in range(1,n):
	seq = filter(lambda subpep:set(get_cyclic_spectrum(subpep)) < set(cyclospec), set(seq))
	if rep != n-1:
		seq = append_c(seq,names)

# format proteins 
cp_seq = ['-'.join([str(int(weights[p])) for p in pep]) for pep in seq]

# print answer and write to out
print(' '.join(cp_seq))
with open('out.txt', 'w') as output:
	output.write(' '.join(cp_seq))