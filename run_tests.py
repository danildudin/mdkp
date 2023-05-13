import sys
import subprocess
from random import randint
from time import time

T_CNT = 1000;

class Solution:
	def __init__(self, path, problem):
		self.problem = problem
		self.z = 0
		self.x = []
		self.read(path)

	def __str__(self):
		return "is_feasible: {0}, z: {1}, x: {2}".format(self.problem.is_feasible(self), self.z, self.x)

	def read(self, path):
		with open(path, "r") as f:
			self.z = int(f.readline())
			f.readline()
			self.x = list(map(int, f.readline().split()))

class Problem:
	def __init__(self, path):
		self.a = []
		self.b = []
		self.c = []
		self.read(path)

	def __str__(self):
		return "c: {0}, b: {1}, a: {2}".format(self.c, self.b, self.a)

	def read(self, path):
		with open(path, "r") as f:
			n, m, real_ans = list(map(int, f.readline().split()))
			self.c = list(map(int, f.readline().split()))
			for i in range(0, m):
				self.a.append(list(map(int, f.readline().split())))
			b = list(map(int, f.readline().split()))

	def is_feasible(self, res):
		for i in range(0, len(self.b)):
			w_sum = 0
			for num in res.x:
				w_sum += self.a[i][num]

			if w_sum > self.b[i]:
				return False
		return True

def run_test(id, path):
	test_in = path + "/{:03d}.in".format(id)
	test_etalon = path + "/{:03d}.out".format(id)
	test_res = path + "/{:03d}.res".format(id)

	print(test_in + "\t", end="")

	cmd_str = "./coral_concurrent.out < {0} > {1}".format(test_in, test_res)
	start_time = time()
	subprocess.run(cmd_str, shell=True)
	print("{:.4f}\t".format(time() - start_time), end="")

	problem = Problem(test_in)
	res_etalon = Solution(test_etalon, problem)
	res = Solution(test_res, problem)

	if res_etalon.z == res.z and problem.is_feasible(res):
		print("ok")
	else:
		print("not ok")
		print("expected: ", res_etalon)
		print("got: ", res)

if __name__ == "__main__":
	path = "tests/small"
	if len(sys.argv) > 1 and sys.argv[1] != "":
		path = sys.argv[1]

	start = 0
	if len(sys.argv) > 2 and sys.argv[2] != "":
		start = int(sys.argv[2])

	end = T_CNT
	if len(sys.argv) > 3 and sys.argv[3] != "":
		end = int(sys.argv[3])

	for i in range(start, end):
		run_test(i, path);
