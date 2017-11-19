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
	chain_info = []
	file = open(input_file)
	parameters = file.readline()
	parameter_list = parameters.split()
	return parameter_list

	
def main():
	parameter_list = import_raw_datas(input_file)
	print(parameter_list)

main()
# if __name__ == '__main__':
# 	main()