import sys
import subprocess
from random import randint

T_CNT = 100;
N_MIN = 10;
N_MAX = 25;
M_MIN = 1;
M_MAX = 8;
C_MIN = 1;
C_MAX = 1000;
W_MIN = 1;
W_MAX = 1000;
B_MIN = 1000;
B_MAX = 6500;

def gen_test(id):
	n = randint(N_MIN, N_MAX)
	m = randint(M_MIN, M_MAX)
	
	c = [randint(C_MIN, C_MAX) for j in range(0, n)]
	b = [randint(B_MIN, B_MAX) for i in range(0, m)]
	a = [[randint(W_MIN, W_MAX) for j in range(0, n)] for i in range(0, m)]

	test_in = "tests/{:02d}.in".format(id)
	test_out = "tests/{:02d}.out".format(id)

	print(test_in)

	with open(test_in, "w") as f:
		f.write("{0} {1}\n".format(n, m))
		for val in c:
			f.write(str(val) + " ")
		f.write("\n")
		for i in range(0, len(a)):
			for val in a[i]:
				f.write(str(val) + " ")
			f.write(str(b[i]))	
			f.write("\n")

	cmd_str = "./mdkp_simple.out < {0} > {1}".format(test_in, test_out)
	subprocess.run(cmd_str, shell=True)


for i in range(0, T_CNT):
	gen_test(i);
