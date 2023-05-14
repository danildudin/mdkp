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

class TestRunner:
	def __init__(self, args):
		self.args = args
		self.exec_time = []
		self.not_ok_cnt = 0

	def run_tests(self):
		regex = None
		if self.args.filter:
			regex = re.compile(self.args.filter)

		for file in sorted(os.listdir(self.args.path)):
			if regex == None or regex.match(file):
				if self.args.list:
					print(file)
				else:
					self.run_test(file)

	def run_test(self, file_name):
		test_in = self.args.path + "/" + file_name
		test_res = test_in + "_res"
		test_etalon = test_in + "_etalon"
		cmd_str = self.args.command + " < {0} > {1}".format(test_in, test_res)

		print(test_in + "\t", end="")

		start_time = time.time()
		subprocess.run(cmd_str, shell=True)
		ts = time.time() - start_time
		self.exec_time.append(ts)

		print("{:.4f}\t".format(ts), end="")

		if self.args.write_etalon:
			subprocess.run("cp {0} {1}".format(test_res, test_etalon), shell=True)
			print()
		elif not self.args.skip_checks:
			problem = Problem(test_in)
			res_etalon = Solution(test_etalon, problem)
			res = Solution(test_res, problem)

			if res_etalon.z == res.z and problem.is_feasible(res):
				print("ok")
			else:
				print("not ok")
				print("expected: ", res_etalon)
				print("got: ", res)
				self.not_ok_cnt += 1
		else:
			print()

	def print_stat(self):
		if args.list:
			return;

		print("tests cnt: {0}".format(len(self.exec_time)))
		print("ok: {0}".format(len(self.exec_time) - self.not_ok_cnt))
		print("not ok: {0}".format(self.not_ok_cnt))
		print("average time (ms): {:04f}".format(sum(self.exec_time) / len(self.exec_time)))
		print("median time (ms): {:04f}".format(sorted(self.exec_time)[len(self.exec_time) // 2]))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Test mdkp.')
	parser.add_argument('-c', '--command', required=True, help='mdkp binary')
	parser.add_argument('-p', '--path', required=True, help='path to the test directory')
	parser.add_argument('-f', '--filter', help='regular expression filter for test files')
	parser.add_argument('-s', '--skip_checks', action='store_true', help='skip comparing the result with the etalon')
	parser.add_argument('-w', '--write_etalon', action='store_true', help='write etalon')
	parser.add_argument('-l', '--list', action='store_true', help='show the list of tests')
	args = parser.parse_args()

	runner = TestRunner(args)
	runner.run_tests()
	runner.print_stat()