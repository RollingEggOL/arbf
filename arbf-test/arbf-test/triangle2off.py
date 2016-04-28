import sys

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print('USAGE: triangle2off <input_node_file> <input ele file> <input edge file> <output OFF file>\n')

	node_file = sys.argv[1]
	ele_file = sys.argv[2]
	edge_file = sys.argv[3]
	off_file = sys.argv[4]

	nv, nf, ne, dim, nodes_per_triangle = 0, 0, 0, 0, 0

	f1 = open(node_file)
	f2 = open(ele_file)
	f3 = open(edge_file)
	f4 = open(off_file, 'w')
	
	# read nv, nf, dim, nodes_per_triangle
	arr = f1.readline().split()
	nv = arr[0]
	dim = arr[1]
	arr = f2.readline().split()
	nf = arr[0]
	nodes_per_triangle = arr[1]
	arr = f3.readline().split()
	ne = arr[0]
	f4.write('OFF\n')
	f4.write('{} {} {}\n'.format(nv, nf, ne))

	for line in f1:
		if not line.startswith('#'):
			arr = line.split()
			f4.write('{} {} 0.0\n'.format(arr[1], arr[2])) # discard node indices
	
	for line in f2:
		if not line.startswith('#'):
			arr = line.split()
			a = int(arr[1]) - 1
			b = int(arr[2]) - 1
			c = int(arr[3]) - 1
			f4.write('{} {} {} {}\n'.format(nodes_per_triangle, a, b, c))
	
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	print('Done!')
