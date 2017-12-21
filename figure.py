#! /usr/bin/env/python
#encoding:utf-8
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

chainsFlie = "newChainInfo.txt"
len_simuBox = 64
maxLen_neigh = 2.5;
# height_simuBox = 64

def builTheBox(len_simuBox,height_simuBox):
	"""
		the function make that building a simulation box 
		and return the API : lst_atomInfo which has record the information of all of atoms;
	"""
	number_totalAtoms = len_simuBox*len_simuBox*height_simuBox
	lst_atomInfo = []
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
				lst_atomInfo.append(singleAtom)			
	if iterNumber != number_totalAtoms:
		print("an error happens in the function of builTheBox")
	print("the total atom number:",number_totalAtoms)
	return lst_atomInfo

def dataImport(chainsFlie): 
	"""
		the function is made to read the datasteam from the original files,
		which has recorded the chains' information about ,chainId and the site coordinates;
		The API
		parameters list and chainOrder list~list;
		but finding the time is to much;
	"""
	list_newChainsInfo = []
	file = open(chainsFlie)
	parameters = file.readline()
	parameterList = parameters.split()
	numOfChain = int(parameterList[-1])
	print("the number of polymer:",numOfChain)
	# the data read from the importfile is the type of "str" need to be translated as a "int" or "float" type.
	iterNumber = 0
	for lines in file.readlines():
		chain_info = lines.split()
		list_newChainsInfo.append(chain_info)  #import the original datas into an list list_newChainsInfo[];
		iterNumber = iterNumber + 1
		pass
	# print(list_newChainsInfo[0])
	if iterNumber != numOfChain:
		print("error in read importFile",iterNumber,numOfChain )
		pass
	file.close()
	print(len(chain_info),list_newChainsInfo[-1])
	return parameterList,list_newChainsInfo

def meanSquareRee2andRg2(parameterList,list_newChainsInfo,lst_atomInfo):
	"""
		this function is used to handle and calculate the properities of 
		the AB diblock copolymer chains,which Sol select A monomer;
		just about the Rg2 and E-E2; 
	"""
	print(parameterList)
	len_Asegment = int(parameterList[3])
	len_Bsegment = int(parameterList[4])
	lenOfPolymer = len_Asegment + len_Bsegment
	# print(len_Asegment,len_Bsegment,lenOfPolymer)
	# print(len(list_newChainsInfo))
	
	# chain = list_newChainsInfo[-1]
	# for x in range(1,len(chain)):
	# 	t = chain[x]
	# 	print(t)
	for chain in list_newChainsInfo:
		nct = 0
		for y in range(1,len(list_newChainsInfo[x])):
			atom_id = list_newChainsInfo[x][y]
			xi = lst_atomInfo[int(atom_id)][1]
			yi = lst_atomInfo[int(atom_id)][2]
			zi = lst_atomInfo[int(atom_id)][3]
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
			lst_atomInfo[int(y)][0] = atom_name
			lst_atomInfo[int(y)][-1] = atom_type
		pass
		if nct != len(list_newChainsInfo[x]):
			print("error in for loop of function meanSquareRee2andRg2")
	nct = nct *len(list_newChainsInfo)	
	print("number of atoms:",nct)
	# check the data correctly;
	num_A =0
	num_B = 0
	for x in range(len(lst_atomInfo)):
		if lst_atomInfo[x][-1] == 190:	
			num_A = num_A + 1
		elif lst_atomInfo[x][-1] == 110:
			num_B = num_B + 1
		pass
	if (num_A+num_B) != nct:
		print("error in function meanSquareRee2andRg2")
	#  calculate the mean squre end-to-end distance
	EE2 = []
	Rg2 = []
	height_simuBox = int(parameterList[1])
	for chain in list_newChainsInfo:
		start_monomer = int(chain[len_Asegment])
		end_monomer = int(chain[-1])
		rx = lst_atomInfo[start_monomer][1]-lst_atomInfo[end_monomer][1]
		ry = lst_atomInfo[start_monomer][2]-lst_atomInfo[end_monomer][2]
		rz = lst_atomInfo[start_monomer][3]-lst_atomInfo[end_monomer][3]
		rz = period_Boundary_condition(rz,height_simuBox)
		rr2 = rx*rx + ry*ry + rz*rz
		EE2.append(rr2)
		pass
	print(EE2)
	meanSquareRee2 = float(sum(EE2)/len(EE2))
	#  program for calculating the meanSquareRg2
	for chain in list_newChainsInfo:
		vector_r2 =[]
		for current_monomer in chain[len_Asegment:]:
			for next_monomer in chain[len_Asegment:]:
				rx = lst_atomInfo[int(current_monomer)][1]-lst_atomInfo[int(next_monomer)][1]
				ry = lst_atomInfo[int(current_monomer)][2]-lst_atomInfo[int(next_monomer)][2]
				rz = lst_atomInfo[int(current_monomer)][3]-lst_atomInfo[int(next_monomer)][3]
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

def build_neighbours_list(lst_atomInfo):
	lst_neighbours = []
	coor  = [];
	near_coor = []
	for x in range(-1,2,1):
		for y in range(-1,2,1):
			for z in range(-1,2,1):
				rr = x**2 + y**2 + z**2
				if rr <maxLen_neigh and rr >0:
					coor = [x,y,z]
					near_coor.append(coor)
					pass
	# print(len(near_coor),maxLen_neigh,near_coor)
	print(lst_atomInfo[0])
	# print(len(lst_atomInfo),64**3)
	for atom in range(0,len(lst_atomInfo)):
		single_neigh_info = []
		for neigh in range(0,len(near_coor)):
			x = lst_atomInfo[atom][1] + near_coor[neigh][0]
			y = lst_atomInfo[atom][2] + near_coor[neigh][1]
			z = lst_atomInfo[atom][3] + near_coor[neigh][2]
			x = period_Boundary_condition(x,len_simuBox,0)
			y = period_Boundary_condition(y,len_simuBox,0)
			z = period_Boundary_condition(z,len_simuBox,0)
			id = len_simuBox**2 *x + len_simuBox*y + z

			single_neigh_info.append([neigh,id])
			pass
		lst_neighbours.append([atom,single_neigh_info])
		pass
	print(lst_neighbours[0])
def period_Boundary_condition(value,period,boundary_judgment):
	"""
		the function of the period_Boundary_condition;
		parameters: value--> the value needs to be judged;
					period--> the period of one direction;
					boundary_judgment --> the judgment for which condition;
						0 >> boundary: like neighbour;
						others >> distance calculating
	"""
	half_period = int(period/2)
	if boundary_judgment == 0:
		boundary_upper_plus = period
		boundary_down_plus = 0
	else :
		boundary_upper_plus = half_period
		boundary_down_plus = -half_period
		pass

	if value > boundary_upper_plus:
		return value - period
	elif value < boundary_down_plus:
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
	parameterList,list_newChainsInfo = dataImport(chainsFlie)
	height_simuBox = int(parameterList[1])
	# print(height_simuBox)
	lst_atomInfo = builTheBox(len_simuBox,height_simuBox)
	# print(lst_atomInfo)
	chosenSegment_type = 2
	# meanSquareRee2andRg2(parameterList,list_newChainsInfo,lst_atomInfo)
	dataExport()
	# build_neighbours_list(lst_atomInfo)
	figures()
	print("main function()")
	pass

main()
# if __name__ == '__main__':
# 	main()