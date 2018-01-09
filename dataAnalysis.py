#! /usr/bin/env/python
#encoding = utf-8
from __future__ import division
import matplotlib.pyplot as plt

importfile = "newChainInfo.txt"

def Box_building(len_simuBox):
	"""
		the function make that building a simulation box
		and return the API : lst_atomInfo which has record the information of all of atoms;
	"""
	number_totalAtoms = len_simuBox**3
	lst_atomInfo = []
	singleAtom = ["Na",0,0,0,0]
	for x in range(0,len_simuBox):
		for y in range(0,len_simuBox):
			for z in range(0,len_simuBox):
				singleAtom = ["Na",x,y,z,0]
				lst_atomInfo.append(singleAtom)			
	return lst_atomInfo

def data_import(rawfile):
	"""
		import data form the "rawfile"
		return the list information of atom,[name,x,y,z,type]
	"""
	lst_parameters = []
	lst_chainInfo = []
	file = open(rawfile)
	parameters_line = file.readline()
	lst_parameters = parameters_line.strip("\n").split(',')
	iterNum = 0
	numOfChain = int(lst_parameters[-3])
	for line in file.readlines():
		line = line.strip('\n')
		lst_chainInfo.append(line.split(','))
		iterNum = iterNum + 1
	if iterNum != numOfChain:
		print("error in read importFile",iterNum,numOfChain )
	return lst_parameters,lst_chainInfo

def dataAnalysis(parameterList,newChainInfo,info_of_atoms):
	len_Asegment = int(parameterList[3])
	len_Bsegment = int(parameterList[4])
	lenOfPolymer = len_Asegment + len_Bsegment
	# print(type(lenOfPolymer),"lenOfPolymer:",lenOfPolymer)
	print(len(newChainInfo))
	for x in range(0,len(newChainInfo)):
		for y in newChainInfo[x][1:-1]:
			xi = info_of_atoms[int(y)][1]
			yi = info_of_atoms[int(y)][2]
			zi = info_of_atoms[int(y)][3]
			ids = xi*len_simuBox*len_simuBox + yi*len_simuBox + zi
			if int(y) != ids:
				# print("error in the transform the format of atoms")
				pass
			pass
		pass

def chainRg2_cal(lst_parameters,lst_atomInfo,lst_chainInfo):
	pass

def chainRee2_cal(lst_parameters,lst_atomInfo,lst_chainInfo):
	for chain in lst_chainInfo:
		dis_end2end = distance(coor1,coor2,len_period)
	pass

def distance(coor1,coor2,len_period):
	'''
		period boundary condition of z direction;\
		coor1 & coor2 is the list of coordinates of two atoms;
		the 3th arguments is the length of period for z direction;
	'''
	x = coor1[0] - coor2[0]
	y = coor1[1] - coor2[1]
	z = coor1[2] - coor2[2]
	# if x > len_period /2:
	# 	x =  len_period - x
	# elif x < -len_period/2:
	# 	x =  len_period + x
	# if y >len_period/2:
	# 	y = len_period -x 
	# elif y < -len_period/2:
	# 	y = len_period + x
	if z > len_period/2:
		z = len_period - z
	elif z < -(len_period/2):
		z = len_period + z
	return x**2 + y **2 + z**2

def plotMayaVi():
	'''
		plot the figure of MayaVi
	'''
	fp = open("B-density.vtk","w")
	fp.close()
	pass

def plotChem3D():
	'''
		plot the figure of ChemBio3D
	'''
	fp = open("B-mono.cc1","w")
	fp.close()
	pass

def main():
	lst_parameters,lst_chainInfo = data_import(importfile)
	lst_atomInfo = Box_building(int(lst_parameters[1]))
	print(lst_parameters)
	print(lst_chainInfo[0])
	plotMayaVi()
	plotChem3D()

main()