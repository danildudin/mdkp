import argparse
import os
import re
import subprocess
import sys
import time

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

def run_test(file_name, args):
	test_in = args.path + "/" + file_name
	test_res = test_in + "_res"
	test_etalon = test_in + "_etalon"
	cmd_str = args.command + " < {0} > {1}".format(test_in, test_res)

	print(test_in + "\t", end="")
	start_time = time.time()
	subprocess.run(cmd_str, shell=True)
	print("{:.4f}\t".format(time.time() - start_time), end="")

	if args.write_etalon:
		subprocess.run("cp {0} {1}".format(test_res, test_etalon), shell=True)
		print()
	elif not args.skip_checks:
		problem = Problem(test_in)
		res_etalon = Solution(test_etalon, problem)
		res = Solution(test_res, problem)

		if res_etalon.z == res.z and problem.is_feasible(res):
			print("ok")
		else:
			print("not ok")
			print("expected: ", res_etalon)
			print("got: ", res)
	else:
		print()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Test mdkp.')
	parser.add_argument('-c', '--command', required=True, help='mdkp binary')
	parser.add_argument('-p', '--path', required=True, help='path to tests directory')
	parser.add_argument('-f', '--filter', help='regexp filter for test files')
	parser.add_argument('-s', '--skip_checks', action='store_true', help='skip comparing result with etalon')
	parser.add_argument('-w', '--write_etalon', action='store_true', help='write answer to etalon')
	parser.add_argument('-l', '--list', action='store_true', help='show list of tests')
	args = parser.parse_args()

	regex = None
	if args.filter:
		regex = re.compile(args.filter)

	for file in os.listdir(args.path):
		if regex == None or regex.match(file):
			if args.list:
				print(file)
			else:
				run_test(file, args)
