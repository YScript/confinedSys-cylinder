#! /usr/bin/env/python
#encoding:utf-8
import numpy as np
import matplotlib.pyplot as plt
import math

chainInfoFile = "newChainInfo.txt";

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
	numOfChain = int(lst_parameters[-2])
	for line in file.readlines():
		line = line.strip('\n')
		lst_chainInfo.append(line.split(','))
		iterNum = iterNum + 1
	if iterNum != numOfChain:
		print("error in read importFile",iterNum,numOfChain )
	return lst_parameters,lst_chainInfo

def buildTheBox(len_simuBox):
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
	# print("the total atom number:",number_totalAtoms)
	return lst_atomInfo

def data_mapping(lst_parameters,lst_atomInfo,lst_chainInfo):
	"""
			the function is used to do with the data to 2D-mapping built with 
		theta*r -- Z.but the procedure only suitable for the one cylce packing;
	"""	
	len_Asegment,len_Bsegment,number_chain,equivalent_radius = lst_parameters[3:]
	len_simuBox = int(lst_parameters[1])
	cylinderradius = float(lst_parameters[0])
	iterNumber = 0
	lst_mapping =[]
	for chain in lst_chainInfo:
		for atom in chain[-int(len_Bsegment):]:
			x,y,z = lst_atomInfo[int(atom)][1:4]
			r = math.sqrt((x- len_simuBox/2)**2+(y- len_simuBox/2)**2)
			if r > cylinderradius:
				print("error in calculating the radius >>> data_mapping()")
				exit(1)
			elif r >float(cylinderradius)/2:   # the value after "r>" need to change manually according different system;
			# else:   # the value after "r>" need to change manually according different system;
				dis_x,dis_y = x- len_simuBox/2,y- len_simuBox/2
				r = math.sqrt((abs(dis_x))**2 + (abs(dis_y))**2)
				if dis_y >= 0:
					theta =math.acos(dis_x/r) 
				else:
					theta = np.pi*2 - math.acos(dis_x/r)
				lst_mapping.append([float(equivalent_radius)*theta,z])
				iterNumber = iterNumber + 1
	if iterNumber != int(len_Bsegment) * int(number_chain):
		print('error in reading the lst_chainInfo || couter the hydrophobic number >>> data_mapping()')
		print("the number of hydrophobic monomer:",iterNumber,int(len_Bsegment)*int(number_chain))
	return lst_mapping

def scatter_figuring(lst_parameters,figureData):
	"""
		scatter_figuring function
	"""
	equivalent_radius = float(lst_parameters[-1])
	xdata = []
	ydata = []
	for data in figureData:
		xdata.append(data[0])
		ydata.append(data[1])
	fig = plt.figure()
	plt.scatter(xdata,ydata,s=30,marker="o")
	plt.title("mapping to (theta*r,Z)")
	plt.axis([0,math.pi*2*equivalent_radius,0,int(lst_parameters[1])])
	plt.xlim({math.pi*2*equivalent_radius,0})
	plt.xticks([0, np.pi*equivalent_radius/2,np.pi*equivalent_radius,3*np.pi*equivalent_radius/2,2*np.pi*equivalent_radius],
		[r'0', r'$\pi r/2$', r'$\pi r$', r'$3\pi r/2$', r'$2\pi r$'])
	plt.xlabel(r"$\theta*r$")
	plt.ylabel("height")
	plt.savefig("mapping.png")
	plt.show()
	pass

def meanSquareRadiusOfGyration(lst_parameters,lst_atomInfo,lst_chainInfo):
	for x in xrange(1,10):
		pass
	pass
def main():
	lst_parameters,lst_chainInfo = data_import(chainInfoFile)
	print("parameters:",lst_parameters)
	lst_atomInfo = buildTheBox(int(lst_parameters[1]))
	lst_mapping = data_mapping(lst_parameters,lst_atomInfo,lst_chainInfo)
	scatter_figuring(lst_parameters,lst_mapping)
	pass

# if __name__ == '__main__':
main()