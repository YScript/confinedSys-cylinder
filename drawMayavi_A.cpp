#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    ifstream fin;
    fin.open("a2.cc1");
    ofstream fout("a2.vtk");

/*    fin.open("step30.cc1");
    ofstream fout("step30.vtk");*/

	int kkk;
	int x, y, z;
    int lx = 60;
    int ly = lx;
    int lz = lx;
    int ntot = lx * ly * lz;
    
    int **atom;
    atom = new int*[3];
    for(int i = 0; i < 3; i++)
    {
            atom[i] = new int[ntot];
    }
    
    int *icha = new int[ntot];
    
    double *dens = new double[ntot];
    
    int **nna;
    nna = new int*[18];
    for(int i = 0; i < 18; i++)
    {
            nna[i] = new int[ntot];
    }
    
    
    
    int offset[18][3] = {1,0,0,
                        -1,0,0,
                        0,1,0,
			0,-1,0,
		0,0,1,
						0,0,-1,
						1,1,0,
						-1,1,0,
						-1,-1,0,
						1,-1,0,
						1,0,1,
						-1,0,-1,
						1,0,-1,
						-1,0,1,
						0,1,1,
						0,-1,-1,
						0,-1,1,
						0,1,-1
                    };
    kkk = 0;                  
    for(int i = 0; i < lx; i++)
    {
    	for(int j = 0; j < ly; j++)
    	{
    		for(int k = 0; k < lz; k++)
    		{
    			atom[0][kkk] = i;
    			atom[1][kkk] = j;
    			atom[2][kkk] = k;
    			icha[kkk] = 0;
    			kkk++;
    		}
    	}
    }
	

    for(int i = 0; i < ntot; i++)
    {
    	for(int j = 0; j < 18; j++)
    	{
    		x = atom[0][i];
    		y = atom[1][i];
    		z = atom[2][i];
    		x += offset[j][0];
    		y += offset[j][1];
    		z += offset[j][2];
    		if(x < 0)
    			x += lx;
    		else if(x >= lx)
    			x = x - lx;
    		else;
    
    		if(y < 0)
    			y += ly;
    		else if(y >= ly)
    			y = y - ly;
    		else;
    		
    		if(z <0)
    			z += lz;
    		else if(z >= lz)
    			z = z - lz;
    		else;
    
    	    nna[j][i] = x*ly*lz + y*lz + z;
    	
    	}
    }
    
 //   cout << "hello world"<<endl; 
 
    
    int n, tmp, c, id;
    double xx, yy, zz;
    string str;
    fin >> n;
  //  cout << n<<endl;
    for(int i = 0; i < n; i++)
    {
            fin >> str;
            fin >> xx;
            fin >> xx;
            fin >> yy;
            fin >> zz;
            fin >> c;
            
            id = (int)xx*ly*lz + (int)yy*lz + (int)zz;
            
            if(c == 110||c == 190)                                  /////////////////////////////////////////////////////////////////
                 icha[id] = 1;
            
            //cout << str << "\t" << xx << "\t" << yy << "\t" << zz << "\t" << c << endl;
            
    }
    fin.close();
//    cout << "hello world"<<endl; 
    for(int i = 0; i < ntot; i++)
    {
            dens[i] = 0.0;
            for(int j = 0; j < 18; j++)
            {
                    if(icha[nna[j][i]] == 1)
                        dens[i] += 1.0;
            }
            dens[i] /= 18;
    }
    



fout << "# vtk DataFile Version 2.0" << endl;
fout << "CT Cylinder interface " << endl;
fout << "ASCII" << endl;
fout << endl;
fout << "DATASET STRUCTURED_POINTS" << endl;
fout << "DIMENSIONS" << "\t" << lx << "\t" << ly << "\t" << lz << endl;
fout << "ORIGIN 0.000000 0.000000 0.000000" << endl;
fout << "SPACING 0.100000 0.100000 0.100000" << endl;
fout << endl;
fout << "POINT_DATA" << "\t" << ntot << endl;
fout << "SCALARS scalars float" << endl;
fout << endl;
fout << "LOOKUP_TABLE default" << endl;
fout << endl;

for(int i = 0; i < ntot; i++)
{
        fout << dens[i] << "\t";
        if(i % (lx*ly) == 0)
             fout << endl;
}


fout.close();
    
    delete []icha;
    for(int i = 0; i < 3; i++)
    {
            delete []atom[i];
    }
    delete []atom;
    
    for(int i = 0; i < 18; i++)
    {
            delete []nna[i];
    }
    delete []nna;
    delete []dens;
    
 //  system("pause");
    return 0;
}
 
