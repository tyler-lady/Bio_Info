# A list Sk+1 of error-free DNA (k+1)-mers (k≤5) taken from the same strand of a circular chromosome (of length ≤50)
# All circular strings assembled by complete cycles in the de Bruijn graph Bk of Sk+1. The strings may be given in any order, but each one should begin with the first (k+1)-mer provided in the input.

# Builds the complete cycle coverings from edges
def cycles(s, edges, k):
	# Determine possible elements to add to the covering
	covering = [i for i, j in enumerate(edges) if j[0] == s[-k+1:]]

	# If nothing left to add, return string if perfect covering (no edges left)
	# If not perfect covering, return empty list
	if len(covering) == 0:
		return [s if edges == [] else []]
	# else move forward with all possible coverings
	else:
		return [cycles(s+edges[x][1][-1], edges[:x]+edges[x+1:], k) for x in covering]

# Generator for unpacking nested lists
def unpack(lst):
	for i in lst:
		if isinstance(i, list):
			for j in unpack(i):
				yield j		
		else:
			yield i

# Main execution of the program
if __name__ == '__main__':
	# Open and read k-mers from our dataset
	with open('dataset.txt') as input_data: #rosalind_grep.txt
		k_mers = [line.strip() for line in input_data.readlines()]

	# Create edges
	k = len(k_mers[0])
	edge = lambda x: [x[0:k-1],x[1:k]]
	edges = [edge(i) for i in k_mers[1:]]

	# Get strings. Flatten and cut the coverings
	# Init coverings with k_mers[0] to start with the first element
	strings = [cycle[:len(k_mers)] for cycle in set(unpack(cycles(k_mers[0], edges, k)))]

	# Print and save the output to new file
	print ('\n'.join(strings))
	with open('097_GREP.txt', 'w') as output:
		output.write('\n'.join(strings))