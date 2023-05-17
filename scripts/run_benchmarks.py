import argparse
import os
import re
import subprocess
import sys
import time

thread_count = [16, 14, 12, 10, 8, 6, 4, 2, 1]

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Test command with a different number of threads')
	parser.add_argument('-c', '--command', required=True)
	parser.add_argument('-p', '--path', required=True)
	args = parser.parse_args()

	for cnt in thread_count: 
		ready_cmd = re.sub(r'--thread_count \d+', '--thread_count {0}'.format(cnt), args.command)
		out_file = args.path + "/" + str(cnt)
		cmd_str = ready_cmd + " > " + out_file

		print(cmd_str + "\t", end="")
		start_time = time.time()
		subprocess.run(cmd_str, shell=True)
		ts = time.time() - start_time
		print("{:.4f}\t".format(ts))
