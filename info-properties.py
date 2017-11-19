#! /usr/bin/env/python3
#encoding:utf-8

from __future__ import division
import numpy as np
import os

input_file = "chain.txt"
len_simuBox = 64

def import_raw_datas(input_file):
	"""
		the function is made to read the datasteam from the original files,
		which has recorded the chains' information about ,chainId and the site coordinates;
		The API
		parameters list and chainOrder list~list;
	"""
	old_chains_info = []
	file = open(input_file)
	parameters = file.readline()
	parameter_list = parameters.split()
	iter_numbers = 0
	for lines in file.readlines():
		chains = lines.split()
		old_chains_info.append(chains)
		iter_numbers = iter_numbers + 1
		pass
	print('the number of chains:',iter_numbers)
	print(old_chains_info)
	return parameter_list

	
def main():
	parameter_list = import_raw_datas(input_file)
	print(parameter_list)

main()
# if __name__ == '__main__':
# 	main()