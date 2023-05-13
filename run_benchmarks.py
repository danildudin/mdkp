import os
import re
import subprocess
import sys
import time

def run_test(test_path):
	test_res = test_path + "_res"

	print(test_path + "\t", end="")
	cmd_str = "./coral_concurrent.out < {0} > {1}".format(test_path, test_res)
	start_time = time.time()
	subprocess.run(cmd_str, shell=True)
	print("{:.4f}\t".format(time.time() - start_time))

if __name__ == "__main__":
	print(sys.argv)

	path = sys.argv[1]

	regex = None
	if len(sys.argv) > 2 and sys.argv[2] != "":
		regex = re.compile(sys.argv[2])

	for file in os.listdir(path):
		if regex == None or regex.match(file):
			run_test(path + "/" + file);
