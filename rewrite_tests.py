import sys

def rewrite_test(id, path):
	test_in = path + "/{:03d}.in".format(id)
	print(test_in)

	b = []
	a = []
	c = []
	with open(test_in, "r") as f:
		n, m = list(map(int, f.readline().split()))
		c = list(map(int, f.readline().split()))
		for i in range(0, m):
			row = list(map(int, f.readline().split()))
			b.append(row[-1])
			row.pop()
			a.append(row)

	with open(test_in, "w") as f:
		f.write("{0} {1} 0\n".format(n, m))
		for val in c:
			f.write(str(val) + " ")
		f.write("\n")
		for i in range(0, len(a)):
			for val in a[i]:
				f.write(str(val) + " ")
			f.write("\n")
		for i in range(0, len(b)):
			f.write(str(b[i]) + " ")
		f.write("\n")

if __name__ == "__main__":
	path = sys.argv[1];
	cnt = int(sys.argv[2])

	for i in range(0, cnt):
		rewrite_test(i, path)
