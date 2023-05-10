import sys
import subprocess
from random import randint

T_CNT = 100;

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
			n, m = list(map(int, f.readline().split()))
			self.c = list(map(int, f.readline().split()))
			for i in range(0, m):
				row = list(map(int, f.readline().split()))
				self.b.append(row[-1])
				row.pop()
				self.a.append(row)

	def is_feasible(self, res):
		for i in range(0, len(self.b)):
			w_sum = 0
			for num in res.x:
				w_sum += self.a[i][num]

			if w_sum > self.b[i]:
				return False
		return True

def run_test(id):
	test_in = "tests/{:02d}.in".format(id)
	test_etalon = "tests/{:02d}.out".format(id)
	test_res = "tests/{:02d}.res".format(id)

	print(test_in + "\t", end="")

	cmd_str = "./mdkp.out < {0} > {1}".format(test_in, test_res)
	subprocess.run(cmd_str, shell=True)

	problem = Problem(test_in)
	res_etalon = Solution(test_etalon, problem)
	res = Solution(test_res, problem)

	if res_etalon.z == res.z and problem.is_feasible(res):
		print("ok")
	else:
		print("not ok")
		print("expected: ", res_etalon)
		print("got: ", res)

for i in range(0, T_CNT):
	run_test(i);
