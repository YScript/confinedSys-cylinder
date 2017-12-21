#include <iostream>
#include <fstream>
/* the c++ produce aim to deal with the rg2 of grafted system;
	respectively , ghp , BCP in the surface micelles , and in the Solutins;	
	made by fengyuan ,2017/4/10.

	using the files "chains.txt";
	
	the data of nchain[][] array is "id_atom -1";
	the produce used in handling the rg2 is repeated, could be pro    moted;
*/
using namespace std;

int const ndx = 3;
int const len_fa = 11;
int const len_fb = 9;
int const len_ghp = 11;
int const number_bcp = 480;
int const zc = 10;
int const hz = 60;
double const BOND2 = 1.6258;

struct coor
{
	int x;
	int y;
	int z;
	int type;
};
int distance_r1_r2(int ,int ,int ,int ,int ,int ,int);

int main(int argc, char const *argv[])
{
	ifstream fin;
	ofstream fout,fo,fp;
	
	fin.open("chains.txt");
    fo.open("chains_surf.txt");
    fp.open("chains_sol.txt");
    fout.open("rg2.txt");

/*	fin.open("step_70.txt");
    fo.open("step_70_surf.txt");
    fp.open("step_70_sol.txt");
    fout.open("step_70_rg2.txt");*/
	int lx = 60, ly = 60 ,lz = hz;	//modify the hz value;
	

	int len_bcp = len_fa + len_fb;
	int lxyz  = lx *ly *lz ;
	int number_ghp = (lx*ly)/(ndx*ndx);
	int number_total_chains = number_ghp + number_bcp;

	int nct = 0;
	coor *atom = new coor[lxyz];
	for (int i = 0; i < lx; ++i)
	{
		for (int j = 0; j < ly; ++j)
		{
			for (int k = 0; k < lz; ++k)
			{
				atom[nct].x = i;
				atom[nct].y = j;
				atom[nct].z = k;
				atom[nct].type = 0;
				nct ++;
			}
		}
	}

	int **nchain;
	nchain = new int*[number_total_chains];
	for (int i = 0; i < number_total_chains; ++i)
	{
		nchain[i] = new int[len_bcp];
	}
	for (int i = 0; i < number_total_chains; ++i)
	{
		for (int j = 0; j < len_bcp; ++j)
		{
			nchain[i][j] = -1;			// easy to test;
		}
	}
	/*the data stream come into three array */
	int *ghp_surf = new int [number_ghp];
	int *bcp_surf = new int [number_bcp];
	int *bcp_sol = new int [number_bcp];

	
	// fin.open("step_30.txt");
	/*string str;
	fin >> str>>ndx;
	fin >> str>>len_fa>>str>>len_fb;
	fin >> str>>number_bcp;

	cout << ndx <<"\t"<<len_fa<<"\t"<<len_fb<<"\t"<<number_bcp<<endl;*/
	int len,id_chain,id_atom;
	for (int i = 0; i < number_total_chains; ++i)
	{
		fin >>id_chain;
		if ( i <number_ghp)
		{
			len = len_ghp;
		}
		else
		{
			len = len_bcp;
		}

		for (int j = 0; j < len; ++j)
		{
			fin >>id_atom;
			nchain[i][j] = id_atom -1;	
		}

	}
	fin.close();
	for (int i = 0; i < number_ghp; ++i)
	{
		ghp_surf[i] = i; 		//record the id of 	ghp_surf;
	}

	bool criterion;	
	int vartest,ncount;
    ncount = 0 ;
    nct = 0;
    for (int i = 0; i < number_bcp; ++i)
    {	
        criterion = true;	
    	id_chain = i + number_ghp;

        for (int j = 0; j < len_fa; ++j){
    	    id_atom = nchain[id_chain][j];
            vartest = atom[id_atom].z;           
            criterion = criterion&&(vartest<=zc && vartest >0);
    	}   	
        if(criterion){
    		bcp_surf[ncount] = id_chain;    //record the id of chain bcp_surf;
    		ncount ++;
    	}else{
            bcp_sol[nct] = id_chain;		//record the id of chain bcp_sol;
            nct ++;
        }
    }
    int number_bcp_surf = ncount;
    int number_bcp_sol = nct;
    
    nct =0;
//    cout <<"number_ghp\t"<<number_ghp<<endl;
//   cout <<"number_bcp_surf\t"<<number_bcp_surf<<endl;
//   cout <<"number_bcp_sol\t"<<number_bcp_sol<<endl;

    double rg2_ghp = 0.0;
    double rg2_bcp = 0.0;
    double rg2_bcp_surf = 0.0;
    double rg2_bcp_sol = 0.0;
    double rgi2;
    int xk,yk,zk,xj,yj,zj;
   
/* ghp_surf*/
    for (int i = 0; i < number_ghp; ++i)
    {
    	rgi2 = 0.0;
    	for (int j = 0; j < len_ghp; ++j)
    	{
    		for (int k = 0; k < len_ghp; ++k)
    		{
    			xj =atom[nchain[i][j]].x;
				yj =atom[nchain[i][j]].y;
				zj =atom[nchain[i][j]].z;

				xk =atom[nchain[i][k]].x;
				yk =atom[nchain[i][k]].y;
				zk =atom[nchain[i][k]].z;

				rgi2 += distance_r1_r2(xj,yj,zj,xk,yk,zk,lx);
    		}
    	}
    	rgi2 /= (double)(2*len_ghp*len_ghp);
    	rg2_ghp +=rgi2; 
    }
    rg2_ghp =(double)(6*rg2_ghp)/(double)(number_ghp*(len_ghp-1)*BOND2);
    cout <<"rg2_ghp\t"<<rg2_ghp<<endl;
/*bcp_surf*/
    for (int i = 0; i < number_bcp_surf; ++i)
    {
    	rgi2 = 0.0;
    	id_chain = bcp_surf[i];
    	for (int j = 0; j < len_fa; ++j)
    	{
    		for (int k = 0; k < len_fa; ++k)
    		{
    			xj =atom[nchain[id_chain][j]].x;
				yj =atom[nchain[id_chain][j]].y;
				zj =atom[nchain[id_chain][j]].z;

				xk =atom[nchain[id_chain][k]].x;
				yk =atom[nchain[id_chain][k]].y;
				zk =atom[nchain[id_chain][k]].z;
				rgi2 += distance_r1_r2(xj,yj,zj,xk,yk,zk,lx);
    		}
    	}
    	rgi2 /= (double)(2*len_fa*len_fa);
    	rg2_bcp_surf += rgi2;
    }
    rg2_bcp_surf =(double)(6*rg2_bcp_surf)/(double)(number_bcp_surf*(len_fa-0.5)*BOND2);
    cout <<"rg2_bcp_surf\t"<<rg2_bcp_surf<<endl;
/* bcp_sol*/
    for (int i = 0; i < number_bcp_sol; ++i)
    {
    	rgi2 = 0.0;
    	id_chain = bcp_sol[i];
    	for (int j = 0; j < len_fa; ++j)
    	{
    		for (int k = 0; k < len_fa; ++k)
    		{
    			xj =atom[nchain[id_chain][j]].x;
				yj =atom[nchain[id_chain][j]].y;
				zj =atom[nchain[id_chain][j]].z;

				xk =atom[nchain[id_chain][k]].x;
				yk =atom[nchain[id_chain][k]].y;
				zk =atom[nchain[id_chain][k]].z;
				rgi2 += distance_r1_r2(xj,yj,zj,xk,yk,zk,lx);
    		}
    	}
    	rgi2 /= (double)(2*len_fa*len_fa);
    	rg2_bcp_sol += rgi2;
    }
    rg2_bcp_sol =(double)(6*rg2_bcp_sol)/(double)(number_bcp_sol*(len_fa-0.5)*BOND2);
    
    rg2_bcp = (rg2_bcp_surf*number_bcp_surf+rg2_bcp_sol*number_bcp_sol)/number_bcp;

   cout <<"rg2_bcp_sol\t"<<rg2_bcp_sol<<endl;
   cout <<"rg2_bcp\t"<<rg2_bcp<<endl;

/*	output the data*/

    fo <<number_bcp_surf<<endl;
    for (int i = 0; i < number_bcp_surf; ++i)
    {
    	id_chain = bcp_surf[i];
    	fo <<id_chain<<"\t";
    	for (int j = 0; j < len_bcp; ++j)
    	{
    		fo <<nchain[id_chain][j]<<"\t";
    	}
    	fo <<endl;
    }
    fp <<number_bcp_sol<<endl;
    for (int i = 0; i < number_bcp_sol; ++i)
    {
    	id_chain = bcp_sol[i];
    	fp <<id_chain<<"\t";
    	for (int j = 0; j < len_bcp; ++j)
    	{
    		fp <<nchain[id_chain][j]<<"\t";
    	}
    	fp <<endl;
    }

    fout <<"rg2_ghp\trg2_bcp_surf\trg2_bcp_sol\trg2_bcp"<<endl;
    fout <<rg2_ghp<<"\t"<<rg2_bcp_surf<<"\t"<<rg2_bcp_sol<<"\t"<<rg2_bcp<<endl;
    if (rg2_bcp_sol >=2.0 || rg2_bcp_surf >= 2.0 || rg2_bcp >=2.0)
    {
    	cout <<"error"<<endl;
    }
    fout.close();
    fo.close();
    fp.close();
/*	free the storage*/
	for (int i = 0; i < number_total_chains; ++i)
	{
		delete []nchain[i];
	}
	delete []nchain;

	delete []atom;
	atom = NULL;
	delete []ghp_surf;
	ghp_surf = NULL;
	delete []bcp_surf;
	bcp_surf = NULL;
	delete []bcp_sol;
	bcp_sol = NULL;
	return 0;
}

int distance_r1_r2(int x1,int y1,int z1,int x2,int y2,int z2,int l){
			int x = x1- x2;
			int y = y1- y2;
			int z = z1- z2;
			int rr2;
			if(x > l/2)
			{	
				x = x -l;
			}
			else if(x < -l/2)
			{
				x = x + l;
			}
			if(y >l/2)
			{
				y = y - l;
			}
			else if(y < -l/2)
			{
				y = y +l;		
			}
			rr2 = x*x + y*y+ z*z;
	return rr2;
}

int Adensity_z_dir(){
	/*
		the subfunction aims to export the data of The Density Distribution of A monomor in Z direction;
		The A monomer identify to class the grafted chain or the free bcp_surf;
	 */
	
}
