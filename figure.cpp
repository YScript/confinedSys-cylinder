#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;
/*
fist,you shoud check the origin input file
	if there need add the code of line 127(fin >>id_of_oldAtom;)
	and line36(fin.getline(buffer,sizeof(buffer)));
	for some old data origin_files;
	subfunction drawMayavi()&rDensityDist()
	
 */
struct coor_property{
	string name;
	int x;
	int y;
	int z;
	int type;
	int color;
};
struct para{
	int number_of_bcp;
	int number_neighs;
	int len_ASegment;
	int len_BSegment;
	double radius;
	double concentration;
	double deta_r;
	int box_lx;
	int box_ly;
	int box_lz;
	int g_lxyz;
} parameters;
double centerRadius;

int main(int argc, char const *argv[]){
// read the parameters from the inputfile;
	int drawMayavi(int number_neighs);
	int radialDensity(struct coor_property *,struct para);
	int mapping(struct coor_property* ,double);
	ifstream fin("chain.txt");
	if (!fin)
	{
		cout <<"error in info;"<<endl;
		exit(1);
	}
	char buffer[256];
	double r0,concentration;
	int cylinderHeight;
	int len_ASegment,len_BSegment;
	int typeA = 1,typeB = 2;
	int colorA = 110,colorB=190;
	int number_of_bcp_import;

	fin >>r0;
	fin >>cylinderHeight;
	fin >>concentration;
	fin >>len_ASegment;
	fin >>len_BSegment;

	int len_of_copolymer = len_ASegment + len_BSegment;
	parameters.len_ASegment = len_ASegment;
	parameters.len_BSegment = len_BSegment;
	int lx = cylinderHeight;	// this needs the var: cylinderHeight read form the import files;
	int ly = lx;
	int lz = lx;
	// build the new space coordinates for mayavi and the old format datas;
	// initialization;
	int g_lxyz = lx*ly*cylinderHeight;
	parameters.box_lx = cylinderHeight;
	parameters.box_ly = cylinderHeight;
	parameters.box_lz = cylinderHeight;
	parameters.g_lxyz = g_lxyz;
	parameters.radius = r0;
	parameters.concentration = concentration;
	coor_property *atom = new coor_property [g_lxyz];

	int counter = 0;
	for (int i = 0; i < lx; ++i){
		for (int j = 0; j < ly; ++j){
			for (int k = 0; k < lz; ++k){
				atom[counter].name = "Ka";
				atom[counter].x = i;
				atom[counter].y = j;
				atom[counter].z = k;
				atom[i].type = 0;
				atom[i].color = 110;
				counter++;
				}
		}
	}
	int x,y,z;
/*build the old space for the initial coordinates and find 
	their identifier of each monomer of each chains;*/

	float dwall = 2.0;
	int cylinderRadius = int(r0+dwall);
	int cylinderDimeters = 2*cylinderRadius+1;
	int number_of_oldParticles = cylinderDimeters*cylinderDimeters*cylinderHeight;
	double rr2;
	coor_property *old_box = new coor_property[number_of_oldParticles];
	counter = 0;
	for (int i = 0; i < cylinderDimeters; ++i){
		for (int j = 0; j < cylinderDimeters; ++j){
			for (int k = 0; k < cylinderHeight; ++k){
				old_box[counter].x = i - cylinderRadius;
				old_box[counter].y = j - cylinderRadius;
				old_box[counter].z = k;
				counter++;
			}
		}
	}
	cout <<"the number of the old box:"<<counter<<endl;
	if (number_of_oldParticles != counter)
	{
		cout <<"error in building the old_box"<<endl;
		cout <<"number_of_oldParticles:\t"<<number_of_oldParticles<<"\t"<<counter<<endl;
		exit(1);
	}
	counter = 0;
	coor_property *old_atom = new coor_property[number_of_oldParticles];
	for (int i = 0; i < number_of_oldParticles; ++i){
		x = old_box[i].x;
		y = old_box[i].y;
		z = old_box[i].z;
		rr2 =sqrt(double(x*x+y*y));
		if (rr2 <= r0){
			old_atom[counter].x = old_box[i].x;
			old_atom[counter].y = old_box[i].y;
			old_atom[counter].z = old_box[i].z;
			counter++;
		}
	}
	int number_of_oldAtom = counter;
	cout <<"the atom number of the system:"<<number_of_oldAtom<<endl;

// read the datas from the inputfile;
	int number_of_bcp = int(number_of_oldAtom*concentration/len_of_copolymer);
	/*if (number_of_bcp != number_of_bcp_import)
	{
		cout <<"error"<<endl;
		exit(1);
	}*/
	int **nchain = new int*[number_of_bcp];
	for (int i = 0; i < number_of_bcp; ++i){
		nchain[i] = new int [len_of_copolymer];
	}

	int id_of_oldAtom = 0;
	int id_of_newAtom = 0;
	int type = 0,color=0;
	string name;
	int id_of_chain;
	counter = 0;
	
	for (int i = 0; i < number_of_bcp; ++i){
		fin >>id_of_chain; // locate this codeline;
		// id_of_chain--;
		for (int j = 0; j < len_of_copolymer; ++j){
			fin >>id_of_oldAtom;
			id_of_oldAtom--;
			nchain[i][j] = id_of_oldAtom;
			if (j <len_ASegment){
				old_atom[nchain[i][j]].name = "Na";
				old_atom[nchain[i][j]].type = typeA;
				old_atom[nchain[i][j]].color = colorA;
			}else{
				old_atom[nchain[i][j]].name = "Br";
				old_atom[nchain[i][j]].type = typeB;
				old_atom[nchain[i][j]].color = colorB;
			}
			x = old_atom[nchain[i][j]].x;
			y = old_atom[nchain[i][j]].y;
			z = old_atom[nchain[i][j]].z;
			id_of_newAtom = int((x+lx/2)*ly*lz+ (y+ly/2)*lz + z);

			nchain[i][j] = id_of_newAtom;
			atom[id_of_newAtom].name = old_atom[id_of_oldAtom].name;
			// atom[id_of_newAtom].x = old_atom[id_of_oldAtom].x;
			// atom[id_of_newAtom].y = old_atom[id_of_oldAtom].y;
			// atom[id_of_newAtom].z = old_atom[id_of_oldAtom].z;
			// the x,y,z coordinates should not change;
			atom[id_of_newAtom].type = old_atom[id_of_oldAtom].type;
			atom[id_of_newAtom].color = old_atom[id_of_oldAtom].color;
			// atom[id_of_newAtom] = old_atom[id_of_oldAtom];
			// counter ++;
		}
	}

	double rr;
	int iterNum = 0;
	/* i have test the number of atoms in the confined System all of the array above
		such as atom[],old_atom[],old_box[]*/
	delete []old_box;
	old_box = NULL;
	delete []old_atom;
	old_atom = NULL;
	fin.close();

// after change the coordinates space and the import nchain info,
	//	we should analyse the datas;

	int number_of_AMonomer = number_of_bcp *len_ASegment;
	int number_of_BMonomer = number_of_bcp *len_BSegment;
	int number_of_solvents;
	cout <<"the number of copolymer:"<<number_of_bcp<<endl;

/* the export partion*/
	ofstream fo,fu,ft,fp;
	fo.open("Amorph.cc1");
	fu.open("Bmorph.cc1");
	ft.open("AB.cc1");
	fp.open("newChainInfo.txt");
	if (!(fo && fu && ft && fp)){
		cout <<"error in open the output file:"<<endl;
		exit(1);
	}

	// dataExport of radial density distribution;
	parameters.deta_r = 0.5;
	radialDensity(atom,parameters);

	// fprintf(A monomer);
	counter = 0;
	fo <<number_of_AMonomer<<endl;
	for (int i = 0; i < g_lxyz; ++i){	
		if (atom[i].type==typeA){
			counter++;
			fo << atom[i].name<<"\t"<<counter<<"\t"<<atom[i].x<<"\t"<<atom[i].y<<"\t"<<atom[i].z<<"\t"<<atom[i].color<<endl;			
		}
	}
	// fprintf(B monomer);
	counter = 0;
	fu <<number_of_BMonomer<<endl;
	for (int i = 0; i < g_lxyz; ++i){
		if (atom[i].type ==typeB){	
			counter ++;
			fu << atom[i].name<<"\t"<<counter<<"\t"<<atom[i].x<<"\t"<<atom[i].y<<"\t"<<atom[i].z<<"\t"<<atom[i].color<<endl;			
		}
	}
	// fprintf(AB, total monomer);
	counter = 0;
	ft <<number_of_AMonomer+number_of_BMonomer<<endl;
	for (int i = 0; i < g_lxyz; ++i){
		if (atom[i].type==typeA || atom[i].type == typeB){
			counter++;
			ft << atom[i].name<<"\t"<<counter<<"\t"<<atom[i].x<<"\t"<<atom[i].y<<"\t"<<atom[i].z<<"\t"<<atom[i].color<<endl;			
		}
	}
	//fprintf(new chains information);
	fp <<r0<<","<<cylinderHeight<<","<<concentration<<","<<len_ASegment<<","<<len_BSegment<<","<<number_of_bcp<<","<<centerRadius<<endl;
	cout <<r0<<"\t"<<cylinderHeight<<"\t"<<concentration<<"\t"<<len_ASegment<<"\t"<<len_BSegment<<"\t"<<number_of_bcp<<endl;
	for (int i = 0; i < number_of_bcp; ++i)
	{
		fp <<i;
		for (int j = 0; j < len_of_copolymer; ++j)
		{
			fp <<","<<nchain[i][j];
		}
		fp <<endl;
	}
	fo.close();
	fu.close();
	ft.close();
	fp.close();
	// close the exportfile;
	// system("pause");

// figure-mayavi-function;
	/*do not creat he neighbour list or array*/
	parameters.number_neighs = 18;
	drawMayavi(parameters.number_neighs);


// free the RAM;
	delete []atom;
	atom = NULL;
	for (int i = 0; i < number_of_bcp; ++i)
	{
		delete []nchain[i];
	}
	delete []nchain;
	nchain = NULL;
	return 0;
}

int drawMayavi(int number_neighs){
	/*
		the subfunction is used to export the .vtk for mayavi-soft;
		needs the neighbour-list to calculate the density of x-y-z; coordinates;
	 */
	cout <<"coming to drawMayavi():"<<endl;
	int **offset = new int*[number_neighs];
	for (int i = 0; i < number_neighs; ++i)
	{
		offset[i] = new int [3];
	}
	// cout <<number_neighs<<endl;
	int counter = 0;
	int x,y,z;
	double rr;
	float len_neighs = 2.5;
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			for (int k = 0; k < 3; ++k){
				x = i-1;
				y = j-1;
				z = k-1;	
				rr = double(x*x+y*y+z*z);
				if ( rr <len_neighs && rr >0 ){
					offset[counter][0] = x;
					offset[counter][1] = y;
					offset[counter][2] = z;
					counter++;
				}
				
			}
		}
	}
	ofstream fout;
	fout.open("Bmorph.vtk");

	// cout <<"hello world"<<endl;

	fout.close();
	return 0;
}

int radialDensity(struct coor_property *atom,struct para parameters){
	/*
		the subfunction is used to export the ".txt" file for the density distribution of 
		radial of the direction
	 */
	// open the ofstream fdenity;
	cout <<"coming to subfunction radialDensity():"<<endl;
	ofstream fdenity,fsample;
	fdenity.open("radialDensity.txt");
	int g_lxyz = parameters.g_lxyz;
	double intervals = parameters.radius /parameters.deta_r ;
	// the intervals control the dimension of the array of radialDensity[] and number_in_unit[];
	double* radialDensity = new double[int(intervals)+1];
	double* number_in_unit = new double[int(intervals)+1];// keep the dimension correctly;
	double x,y,rr;
	int n;
	int iterNum = 0;
	int counter = 0;
	for (int i = 0; i <= intervals; ++i){
		radialDensity[i] = 0.0;
		number_in_unit[i] = 0.0;
	}
	for (int i = 0; i < g_lxyz; ++i)
	{
		x = atom[i].x - (parameters.box_lx/2);
		y = atom[i].y - (parameters.box_ly/2);
		rr = sqrt(x*x + y*y);
		if (rr <= (parameters.radius))
		{
			iterNum ++;
		}
		if (rr == 0.0)
		{
			counter ++;
		}
	}
	cout <<"the atom number in confined Cylinder:"<<iterNum<<"\t"<<endl;
	cout <<"radius:\t"<<parameters.radius<<"\t"<<parameters.deta_r<<"\t"<<intervals<<endl;

	iterNum = 0;
	counter = 0;
	for (int j = 0; j < int(intervals); ++j){
		for (int i = 0; i < g_lxyz; ++i){
			x = atom[i].x - (parameters.box_lx/2);
			y = atom[i].y - (parameters.box_ly/2);
			rr = x*x + y*y;
			if (  sqrt(rr) > j*parameters.deta_r && sqrt(rr) <=(j+1)*parameters.deta_r ){
				number_in_unit[j] += 1.0;
				if (atom[i].type == 2)
				{
					radialDensity[j] += 1.0;
					iterNum ++;
				}
					counter ++;
			}
			if (j== 0 && sqrt(rr) >= 0.0 && sqrt(rr) < parameters.deta_r){
				number_in_unit[0] += 1.0;
				if (atom[i].type == 2)
				{
					radialDensity[0] += 1.0;
					iterNum ++;	
				}
				counter ++;
			}
		}
	}
	cout <<"number_of_BMonomer:"<<iterNum<<endl;
	int number_of_BMonomer = iterNum;
	cout <<"number_of monomer:"<<counter<<endl;

	cout <<"iterNum:\t"<<iterNum<<endl;
	iterNum = 0;
	double newcounter = 0.0;
	double centerPosition=0.0;
	for (int i = 0; i < int(intervals); ++i)
	{
		radialDensity[i] /=number_in_unit[i];
		if (radialDensity[i] > 0)
		{
			centerPosition +=(i+1)*parameters.deta_r;
			iterNum ++;
		}
	}
	centerRadius = centerPosition/iterNum;
	iterNum = 0;
	for (int i = 0; i < int(intervals); ++i)
	{
		iterNum += int(number_in_unit[i]);
		newcounter += radialDensity[i]*number_in_unit[i];
		// cout <<(i+1)*parameters.deta_r <<"\t"<<number_in_unit[i]<<"\t"<<radialDensity[i]<<endl;
		fdenity <<(i+1)*parameters.deta_r <<"\t"<<number_in_unit[i]<<"\t"<<radialDensity[i]<<endl;
	}
	if (newcounter != number_of_BMonomer)
	{
		cout <<"error in counter B monomer in subfunction-radialDensity():"<<counter<<"\t"<<number_of_BMonomer<<endl;
		exit(1);
	}
	cout <<"count the number of volume:"<<iterNum<<endl;
	
	// close the file export stream;
	fdenity.close();
	cout <<"Finishing function radialDensity()"<<endl;
	delete []radialDensity;
	radialDensity = NULL;
	delete []number_in_unit;
	number_in_unit = NULL;
}
