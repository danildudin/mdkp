from random import randint
import argparse
import os
import subprocess
import sys

def gen_test(id, args):
	n = randint(args.n_min, args.n_max)
	m = randint(args.m_min, args.m_max)

	c = [randint(args.c_min, args.c_max) for j in range(0, n)]
	b = [randint(args.b_min, args.b_max) for i in range(0, m)]
	a = [[randint(args.w_min, min(b[i], args.w_max)) for j in range(0, n)] for i in range(0, m)]

	test_in = args.path + "/{:04d}".format(id)
	test_etalon = test_in + "_etalon"

	print(test_in)

	with open(test_in, "w") as f:
		f.write("{0} {1} 0\n".format(n, m))
		for val in c:
			f.write(str(val) + " ")
		f.write("\n")
		for i in range(0, len(a)):
			for val in a[i]:
				f.write(str(val) + " ")
			f.write("\n")
		for i in range(0, m):
			f.write(str(b[i]) + " ")

	cmd_str = args.command + "< {0} > {1}".format(test_in, test_etalon)
	subprocess.run(cmd_str, shell=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Generate tests.')
	parser.add_argument('-c', '--command', required=True, help='mdkp binary')
	parser.add_argument('-p', '--path', required=True, help='path to the test directory')
	parser.add_argument('--tests_cnt', type=int, required=True, help='number of tests')
	parser.add_argument('--n_min', type=int, default=10)
	parser.add_argument('--n_max', type=int, default=10)
	parser.add_argument('--m_min', type=int, default=2)
	parser.add_argument('--m_max', type=int, default=2)
	parser.add_argument('--c_min', type=int, default=1)
	parser.add_argument('--c_max', type=int, default=100)
	parser.add_argument('--w_min', type=int, default=1)
	parser.add_argument('--w_max', type=int, default=100)
	parser.add_argument('--b_min', type=int, default=1)
	parser.add_argument('--b_max', type=int, default=600)
	args = parser.parse_args()

	os.makedirs(args.path, exist_ok=True)

	for i in range(0, args.tests_cnt):
		gen_test(i, args);
