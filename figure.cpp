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
 */
struct coor_property{
	string name;
	int x;
	int y;
	int z;
	int type;
	int color;
};

int main(int argc, char const *argv[]){
	// read the parameters from the inputfile;
	int drawMayavi(int number_neighs);
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

	cout <<"line"<<endl;
	
	fin >>r0;
	fin >>cylinderHeight;
	fin >>concentration;
	fin >>len_ASegment;
	fin >>len_BSegment;
	cout <<r0<<"\t"<<cylinderHeight<<"\t"<<concentration<<"\t"<<len_ASegment<<"\t"<<len_BSegment<<endl;

	int len_of_copolymer = len_ASegment + len_BSegment;
	int lx = cylinderHeight;
	
	int ly = lx;
	int lz = lx;
	// build the new space coordinates for mayavi and the old format datas;
	int g_lxyz = lx*ly*cylinderHeight;
	coor_property *atom = new coor_property [g_lxyz];
	for (int i = 0; i < g_lxyz; ++i){
		atom[i].type = 0;
		atom[i].color = 110;
	}
	int counter = 0;
	for (int i = 0; i < lx; ++i){
		for (int j = 0; j < ly; ++j){
			for (int k = 0; k < lz; ++k){
				atom[counter].name = "Ka";
				atom[counter].x = i;
				atom[counter].y = j;
				atom[counter].z = k;
				counter++;
			}
		}
	}
	/*build the old space for the initial coordinates and 
	 find their identifier of each monomer of evry chains;*/
	float dwall = 2.0;
	int cylinderRadius = int(r0+dwall);
	int cylinderDimeters = 2*cylinderRadius+1;
	int number_of_oldParticles = cylinderDimeters*cylinderDimeters*cylinderHeight;
	int x,y,z;
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
	number_of_oldParticles = counter;
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
	cout <<"the number of the old system:"<<number_of_oldAtom<<endl;

	// read the datas from the inputfile;
	int number_of_bcp = int(number_of_oldAtom*concentration/len_of_copolymer);
	int **nchain = new int*[number_of_bcp];
	for (int i = 0; i < number_of_bcp; ++i){
		nchain[i] = new int [len_of_copolymer];
	}
	
	int id_of_oldAtom = 0;
	int id_of_newAtom = 0;
	int type = 0,color=0;
	string name;
	int id_of_chain;
	for (int i = 0; i < number_of_bcp; ++i){
		fin >>id_of_chain; // locate this codeline;
		for (int j = 0; j < len_of_copolymer; ++j){
			fin >>id_of_oldAtom;
			nchain[i][j] = id_of_oldAtom -1;
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

			atom[id_of_newAtom].name = old_atom[nchain[i][j]].name;
			atom[id_of_newAtom].x = old_atom[nchain[i][j]].x;
			atom[id_of_newAtom].y = old_atom[nchain[i][j]].y;
			atom[id_of_newAtom].z = old_atom[nchain[i][j]].z;
			atom[id_of_newAtom].type = old_atom[nchain[i][j]].type;
			atom[id_of_newAtom].color = old_atom[nchain[i][j]].color;
		}
	}
	delete []old_box;
	old_box = NULL;
	delete []old_atom;
	old_atom = NULL;
	fin.close();

	// after change the coordinates space, we should analyse the datas;
	int number_of_AMonomer = number_of_bcp *len_ASegment;
	int number_of_BMonomer = number_of_bcp *len_BSegment;
	int number_of_solvents;
	cout <<"the number of copolymer:"<<number_of_bcp<<endl;
	cout <<"the number_of_AMonomer:"<<number_of_AMonomer<<endl;
	cout <<"the number_of_BMonomer:"<<number_of_BMonomer<<endl;
	ofstream fo,fu,ft,fp;
	fo.open("Amorph.cc1");
	fu.open("Bmorph.cc1");
	ft.open("AB.cc1");
	if (!(fo && fu && ft && fp)){
		cout <<"error in open the output file:"<<endl;
		exit(1);
	}
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
	fo.close();
	fu.close();
	ft.close();
	fp.close();
	int number_neighs = 18;
	cout <<"coming to drawMayavi()"<<endl;
	drawMayavi(number_neighs);
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
	int **offset = new int*[number_neighs];
	for (int i = 0; i < number_neighs; ++i)
	{
		offset[i] = new int [3];
	}
	cout <<number_neighs<<endl;
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

	cout <<"hello world"<<endl;

	fout.close();
	return 0;
}
