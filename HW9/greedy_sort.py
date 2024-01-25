from operator import neg

def g_sort(p):
    p_seq = []
    # Lambda to find index
    _index = lambda x, k: list(map(abs, x)).index(k)
    # Lambda to reverse sequence
    _sort = lambda x, i, j: x[:i] + list(map(neg, x[i:j+1][::-1])) + x[j+1:]
    p_list = list(p)
    # Loop over permutation to sort
    i = 0
    while i < len(p_list):
        if p_list[i] == i+1:
            i += 1
        else:
            p_list = _sort(p_list, i, _index(p_list, i+1))
            p_seq.append(p_list)
    return p_seq

if __name__ == '__main__':
    # Read input
    with open('dataset.txt') as input_data:
        # Create map of data (later converted to list for manipulation)
        p = map(int, input_data.read().strip().lstrip('(').rstrip(')').split())
    # Get list of reversals
    reversals = g_sort(p)
    # Format the permutations
    reversals = ['('+' '.join([['', '+'][v > 0] + str(v) for v in p])+')' for p in reversals]
    # Write to output
    print('\n'.join(reversals))
    with open('out.txt', 'w') as output_data:
        output_data.write('\n'.join(reversals))