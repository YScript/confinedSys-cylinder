#! /usr/bin/env/python
#encoding:utf-8
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

chainsFlie = "newChainInfo.txt"
info_raw_chain = "chain.txt"
len_simuBox = 64
# height_simuBox = 64

def builTheBox(len_simuBox,height_simuBox):
	"""
		the function make that building a simulation box 
		and return the API : info_of_atoms which has record the information of all of atoms;
	"""
	number_totalAtoms = len_simuBox*len_simuBox*height_simuBox
	info_of_atoms = []
	iterNumber = 0
	singleAtom = ["Sol",0,0,0,0]
	for x in range(0,len_simuBox):
		for y in range(0,len_simuBox):
			for z in range(0,height_simuBox):
				# singleAtom[0] = "Sol"  ## the comments code does not useful;
				# singleAtom[1] = x 	## because of it reference delivery of list in python;
				# singleAtom[2] = y
				# singleAtom[3] = z
				# singleAtom[4] = 0
				singleAtom = ["Sol",x,y,z,0]
				iterNumber = iterNumber + 1
				info_of_atoms.append(singleAtom)			
	if iterNumber != number_totalAtoms:
		print("an error happens in the function of builTheBox")
	print("the total atom number:",number_totalAtoms)
	return info_of_atoms

def dataImport(chainsFlie): 
	"""
		the function is made to read the datasteam from the original files,
		which has recorded the chains' information about ,chainId and the site coordinates;
		The API
		parameters list and chainOrder list~list;
	"""
	newChainInfo = []
	file = open(chainsFlie)
	parameters = file.readline()
	parameterList = parameters.split()
	# confinedRadius = float(parameterList[0])
	# confinedHeight = int(parameterList[1])
	# polymerConcentration = float(parameterList[2])
	# Asegment = int(parameterList[3])
	# Bsegment = int(parameterList[4])
	# lenOfPolymer = Asegment + Bsegment
	# print(confinedRadius,confinedHeight,polymerConcentration,Asegment,Bsegment)
	datas = file.readline()
	dataList = datas.split()
	numOfChain = int(dataList[-1])
	parameterList.append(numOfChain)
	print("the number of polymer:",numOfChain)
	print(parameterList)
	# the data read from the importfile is the type of "str" need to be translated as a "int" or "float" type.
	iterNumber = 0
	for lines in file.readlines():
		chain_info = lines.split()
		newChainInfo.append(chain_info)  #import the original datas into an list newChainInfo[];
		iterNumber = iterNumber + 1
		pass
	if iterNumber != numOfChain:
		print("error in read importFile" )
		pass
	file.close()
	print(len(chain_info),newChainInfo[-1])
	return parameterList,newChainInfo

def read_chainInfo_from_file(info_raw_chain):
	"""
		read the chains information from the raw chain producted by the .f90 procedure;
		the info_raw_chain('chain.txt') does not have the format of 'one chain -- one line '
		so it need to be import as one by one;
	"""
	old_chain_info = []
	file = open(info_raw_chain)
	parameters = file.readline()
	parameterList = parameters.split()
	# the parameterList: r0,height,concentration,len-ASegment,len-BSegment;
	len_Asegment = int(parameterList[4])
	len_Bsegment = int(parameterList[5])
	len_copolymer = len_Asegment + len_Bsegment

	number_totalChains = int(parameterList[3])
	for chain in range(number_totalChains):
		for monomer in range(len_copolymer):
			if monomer < len_Asegment:
				pass
			elif monomer < len_copolymer:
				pass
			else:
				exit(1)
	pass

def meanSquareRee2andRg2(parameterList,newChainInfo,info_of_atoms):
	"""
		this function is used to handle and calculate the properities of 
		the AB diblock copolymer chains,which Sol select A monomer;
		just about the Rg2 and E-E2; 
	"""
	len_Asegment = int(parameterList[3])
	len_Bsegment = int(parameterList[4])
	lenOfPolymer = len_Asegment + len_Bsegment
	# print(len_Asegment,len_Bsegment,lenOfPolymer)
	# print(len(newChainInfo))
	for x in range(0,len(newChainInfo)):
		nct = 0
		for y in newChainInfo[x]:
			xi = info_of_atoms[int(y)][1]
			yi = info_of_atoms[int(y)][2]
			zi = info_of_atoms[int(y)][3]
			ids = xi*len_simuBox*len_simuBox + yi*len_simuBox + zi
			# conform the data correctly;
			if int(y) != ids:
				print("error in the transform the format of atoms")
				pass
			# judge the atom is A or B monomer;
			if nct < (len_Asegment):
				atom_name = "A"
				atom_type= 190
			else:
				atom_name = "B"
				atom_type = 110
			nct = nct + 1
			info_of_atoms[int(y)][0] = atom_name
			info_of_atoms[int(y)][-1] = atom_type
		pass
		if nct != len(newChainInfo[x]):
			print("error in for loop of function meanSquareRee2andRg2")
	nct = nct *len(newChainInfo)	
	print("number of atoms:",nct)
	# check the data correctly;
	num_A =0
	num_B = 0
	for x in range(len(info_of_atoms)):
		if info_of_atoms[x][-1] == 190:	
			num_A = num_A + 1
		elif info_of_atoms[x][-1] == 110:
			num_B = num_B + 1
		pass
	if (num_A+num_B) != nct:
		print("error in function meanSquareRee2andRg2")
	#  calculate the mean squre end-to-end distance
	EE2 = []
	Rg2 = []
	height_simuBox = int(parameterList[1])
	for chain in newChainInfo:
		start_monomer = int(chain[len_Asegment])
		end_monomer = int(chain[-1])
		rx = info_of_atoms[start_monomer][1]-info_of_atoms[end_monomer][1]
		ry = info_of_atoms[start_monomer][2]-info_of_atoms[end_monomer][2]
		rz = info_of_atoms[start_monomer][3]-info_of_atoms[end_monomer][3]
		rz = period_Boundary_condition(rz,height_simuBox)
		rr2 = rx*rx + ry*ry + rz*rz
		EE2.append(rr2)
		pass
	print(EE2)
	meanSquareRee2 = float(sum(EE2)/len(EE2))
	#  program for calculating the meanSquareRg2
	for chain in newChainInfo:
		vector_r2 =[]
		for current_monomer in chain[len_Asegment:]:
			for next_monomer in chain[len_Asegment:]:
				rx = info_of_atoms[int(current_monomer)][1]-info_of_atoms[int(next_monomer)][1]
				ry = info_of_atoms[int(current_monomer)][2]-info_of_atoms[int(next_monomer)][2]
				rz = info_of_atoms[int(current_monomer)][3]-info_of_atoms[int(next_monomer)][3]
				rz = period_Boundary_condition(rz,height_simuBox)
				rr2 = rx*rx + ry*ry + rz*rz
				vector_r2.append(rr2)
				pass
			pass
		Rg2_of_singleChain = float(sum(vector_r2)/(2*len(chain)*len(chain)))
		Rg2.append(Rg2_of_singleChain)
		pass
	meanSquareRg2 = float(sum(Rg2)/len(Rg2))
	if len(Rg2) != len(EE2):
		print("error in calculating Rg2&Ree2")
	print("meanSquareRee2:",meanSquareRee2)
	print("meanSquareRg2:",meanSquareRg2)

def convert2CylidericalCoor():
	pass


def period_Boundary_condition(value,period):
	"""
		the function of the period_Boundary_condition;
		parameters: value--> the value needs to be judged;
					period--> the period of one direction;
	"""
	half_period = int(period/2)
	if value >half_period:
		return value - period
	elif value <= -half_period:
		return value + period
	else:
		return value
def dataExport():
	rg2 = []
	rg2_sum = 0
	# open()
	print("dataExport()")
	pass

def figures():
	print("figures()")
	pass

def main():
	"""
		the main program
	"""
	parameterList,newChainInfo = dataImport(chainsFlie)
	height_simuBox = int(parameterList[1])
	# print(height_simuBox)
	info_of_atoms = builTheBox(len_simuBox,height_simuBox)
	# print(info_of_atoms)
	chosenSegment_type = 2
	meanSquareRee2andRg2(parameterList,newChainInfo,info_of_atoms)
	dataExport()
	figures()
	print("main function()")
	pass

main()
# if __name__ == '__main__':
# 	main()