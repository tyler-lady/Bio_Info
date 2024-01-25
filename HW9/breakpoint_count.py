def get_breakpoints(p):
    p_list = list(p)
    return sum(map(lambda i,j: i - j != 1, p_list+[len(p_list)+1], [0]+p_list))

if __name__ == '__main__':
    # Read input
    with open('dataset.txt') as input:
        p = map(int, input.read().strip().lstrip('(').rstrip(')').split())
    # Get breakpoint count
    breakpoints = get_breakpoints(p)
    # Handle output
    print(str(breakpoints))
    with open('out.txt', 'w') as output_data:
        output_data.write(str(breakpoints))