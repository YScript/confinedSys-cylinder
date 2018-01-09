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
		lst_parameters=[radius_cylinder,height_cylinder,concentration,len_Asegment,len_Bsegment,number_chain,center_position]
	"""
	lst_parameters = []
	lst_chainInfo = []
	file = open(rawfile,'r')
	parameters_line = file.readline() #parameters line;
	lst_parameters = list(map(eval,parameters_line.strip("\n").split(',')))
	iterNum = 0
	for line in file.readlines():
		line = list(map(int,line.strip('\n').split(','))) # map()function return a struct of type "map";
		lst_chainInfo.append(line)
		iterNum = iterNum + 1
	file.close()
	if iterNum != lst_parameters[-2]:
		print("error in read importFile",iterNum,lst_parameters[-2])
	print("parameters:",tuple(lst_parameters))
	print('lst_chainInfo[0]:',lst_chainInfo[0])
	return tuple(lst_parameters),tuple(lst_chainInfo) # the lst_parameters and the lst_chainInfo can't change

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
		for atom in chain[-len_Bsegment:]:
			x,y,z = lst_atomInfo[atom][1:4]
			r = math.sqrt((x- len_simuBox/2)**2+(y- len_simuBox/2)**2)
			if r > cylinderradius:
				print("error in calculating the radius >>> data_mapping()")
				exit(1)
			elif r > 0:
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
	plt.axis([0,math.pi*2*equivalent_radius,0,lst_parameters[1]])
	plt.xlim({0,math.pi*2*equivalent_radius})
	plt.xticks([0, np.pi*equivalent_radius/2,np.pi*equivalent_radius,3*np.pi*equivalent_radius/2,2*np.pi*equivalent_radius],
		[r'0', r'$\pi r/2$', r'$\pi r$', r'$3\pi r/2$', r'$2\pi r$'])
	plt.xlabel(r"$\theta*r$")
	plt.ylabel("height")
	plt.savefig("mapping.png")
	# plt.show()
	pass

def Rg2_calculating(monomer,lst_parameters,lst_atomInfo,lst_chainInfo):
	'''
		calculating the mean square Radius of Gyration of polymer chains;
		lst_parameters=[radius_cylinder,height_cylinder,concentration,len_Asegment,len_Bsegment,center_position]
	'''
	len_Asegment,len_Bsegment = lst_parameters[3:5]
	len_copolymer = len_Asegment + len_Bsegment
	number_chain = lst_parameters[-2]
	height_cylinder = lst_parameters[1]
	
	startPosition = -len_Bsegment
	print(len_Bsegment)
	endPosition = -1
	sum_dist2_gyration = 0.0
	lst_rg2 = []
	Rg2 = 0
	Rg2_ideal = 0
	Rg2_reduce = 0
	for chain in lst_chainInfo:
		x,y,z = 0,0,0
		for atom in chain[startPosition:]:
			# massCenter = list(map(lambda x:x[0]+x[1],zip(massCenter,lst_atomInfo[atom][1:4])))
			x = x + lst_atomInfo[atom][1]
			y = y + lst_atomInfo[atom][2]
			z = z + lst_atomInfo[atom][3]
		print(x,y,z)
		massCenter = list(map(lambda x:x/len_Bsegment , [x,y,z]))#[x,y,z]/len_Bsegment
			# massCenter = list(map(lambda x: x/len_Bsegment,massCenter))
		print(massCenter,len(massCenter))
		sum_dist2 = 0
		for atom in chain[startPosition:]:
			sum_dist2 = sum_dist2 + distanceSquare(lst_atomInfo[atom][1:4],massCenter,height_cylinder)
		
		sum_dist2 = sum_dist2/len_Bsegment
		lst_rg2.append([sum_dist2,massCenter])

	# print(lst_rg2)
	print(massCenter)
	return Rg2,Rg2_ideal,Rg2_reduce

def Ree2_calculating(monomer,lst_parameters,lst_atomInfo,lst_chainInfo):
	'''
		calculating the distance of end-to-end of polymer Chains;
		lst_parameters=[radius_cylinder,height_cylinder,concentration,len_Asegment,len_Bsegment,number_chain,center_position]
	'''
	len_Asegment,len_Bsegment = lst_parameters[3:5]
	len_copolymer = len_Asegment + len_Bsegment
	number_chain = lst_parameters[-2]
	height_cylinder = lst_parameters[1]
	dist_end2end = 0.0
	if monomer == 'A':
		startPosition = 1
		endPosition = len_Asegment
	elif monomer == 'B':
		startPosition = -len_Bsegment
		endPosition = -1
	else:
		startPosition = 1
		endPosition = -1
	sum_dist2_end2end = 0.0
	for chain in lst_chainInfo:
		dist_end2end = distanceSquare(lst_atomInfo[chain[startPosition]][1:4],\
			lst_atomInfo[chain[endPosition]][1:4],height_cylinder)
		sum_dist2_end2end = sum_dist2_end2end + dist_end2end
	Ree2_cal = sum_dist2_end2end/number_chain 
	# print(type(chain),Rstart,Rend,lst_atomInfo[Rstart],lst_atomInfo[Rend])
	print("ree2-"+monomer+':',Ree2_cal)
	return Ree2_cal

def distanceSquare(lst_r1,lst_r2,period_length):
	x = lst_r1[0] - lst_r2[0]
	y = lst_r1[1] - lst_r2[1]
	z = lst_r1[2] - lst_r2[2]
	if z >period_length/2:
		z = -z + period_length
	elif z <-period_length/2:
		z = z + period_length
	r2 = x**2 +y**2 + z**2
	return r2

def main():
	tuple_parameters,tuple_chainInfo = data_import(chainInfoFile)
	lst_atomInfo = buildTheBox(tuple_parameters[1])
	print("lst_atomInfo[0]:",lst_atomInfo[0])
	print(tuple_parameters[4])
	lst_mapping = data_mapping(tuple_parameters,lst_atomInfo,tuple_chainInfo)
	scatter_figuring(tuple_parameters,lst_mapping)
	Ree2_B = Ree2_calculating('B',tuple_parameters,lst_atomInfo,tuple_chainInfo)
	Rg2,Rg2_ideal,Rg2_reduce = Rg2_calculating('B',tuple_parameters,lst_atomInfo,tuple_chainInfo) 
	print("done")
# if __name__ == '__main__':
main()