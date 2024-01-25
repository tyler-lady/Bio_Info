# returns a dict for protein to mass
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

def get_alt_weight_dict():
	table ='''A   71
	C   103
	D   115
	E   129
	F   147
	G   57
	H   137
	I   113
	K   128
	L   113
	M   131
	N   114
	P   97
	Q   128
	R   156
	S   87
	T   101
	V   99
	W   186
	Y   163''' 

	weight_dict = dict()

	for p in table.split('\n'):
		weight_dict[p.strip('\t').split()[0]] = float(p.strip('\t').split()[1])

	return weight_dict

def get_cyclic_spectrum(pep):
	# Dict for RNA to protein
	weight = get_alt_weight_dict()

	# Init as mass 0 and mass of peptide
	spec = [0, sum([int(weight[p]) for p in pep])]

	# Find masses of adjacent subpeptides
	spec += [sum([int(weight[p]) for p in (pep*2)[j:j+i]]) for i in range(1,len(pep)) for j in range(len(pep))]

	# Sort list ascending order and change to strings
	spec = map(str,sorted(spec))

	return spec

if __name__ == '__main__':
	with open('dataset.txt') as input_data:
		pep = input_data.read().strip()

	spec = get_cyclic_spectrum(pep)

	# Print and save answer
	print(' '.join(spec))
	with open('out.txt', 'w') as output_data:
		output_data.write(' '.join(spec))