
#include <cstddef>

namespace nglib{
#include "nglib.h"
}

#include <cstdio>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace nglib;

int main(int argc, char **argv)
{
	string infile="object.stl";
	char file[infile.length()+1];
	strcpy(file,infile.c_str());

	int np, ne, se;
	Ng_Init();
	Ng_Mesh *mesh = Ng_NewMesh();
	Ng_STL_Geometry *stl_geom = Ng_STL_LoadGeometry(file);
	if (!stl_geom) {
		std::cout << "Error reading in STL File:\n";
		return 1;
	}
	std::cout << "Successfully loaded STL File\n";

	Ng_Meshing_Parameters mp;
	mp.maxh = 1.0e+6;
	mp.fineness = 0.4;
	mp.second_order = 0;

	std::cout << "Initialise the STL Geometry structure....\n";
	Ng_Result ng_res = Ng_STL_InitSTLGeometry(stl_geom);
	if (ng_res != NG_OK) {
		std::cout << "Error Initialising the STL Geometry....Aborting!!\n";
		return 1;
	}

	std::cout << "Start Edge Meshing....\n";
	ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
	if (ng_res != NG_OK) {
		std::cout << "Error in Edge Meshing....Aborting!!\n";
		return 1;
	}

	std::cout << "Start Surface Meshing....\n";
	ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
	if (ng_res != NG_OK) {
		std::cout << "Error in Surface Meshing....Aborting!!\n";
		return 1;
	}

	std::cout << "Start Volume Meshing....\n";
	ng_res = Ng_GenerateVolumeMesh(mesh, &mp);
	if (ng_res != NG_OK) {
		std::cout << "Error in Volume Meshing....Aborting!!\n";
		return 1;
	}

	std::cout << "Meshing successfully completed....!!\n";

	// volume mesh output
	np = Ng_GetNP(mesh);
	std::cout << "Points: " << np << "\n";

	ne = Ng_GetNE(mesh);
	std::cout << "Elements: " << ne << "\n";

	std::cout << "Saving Mesh in VOL Format....\n";
	Ng_SaveMesh(mesh, "test.vol");
	se=Ng_GetNSE(mesh);
	// vtk file
	string soubor="out.vtk";
	ofstream sfile;
	sfile.open(soubor);
	sfile << "# vtk DataFile Version 2.0"<< endl;
	sfile <<"NGsolve, Created by NGsolve"<<endl;
	sfile <<"ASCII"<<endl;
	sfile <<"DATASET UNSTRUCTURED_GRID"<<endl;
	sfile <<"POINTS "<<np<<" double"<<endl;
	for (int i=1;i<=np;i++){
		double *point= new double[3];
		Ng_GetPoint(mesh,i,point);
			   sfile <<setprecision(20)<<point[0]<< " "<<point[1] << " "<< point[2]<<endl;
		}
	sfile <<endl;
	Ng_Surface_Element_Type a[se];
	int h=0;
	for (int i=1;i<=se;i++){
		int *element = new int[8];
		a[i-1] = Ng_GetSurfaceElement(mesh,i,element);
		if(a[i-1]==1)h=h+3;
		if(a[i-1]==2)h=h+4;
		if(a[i-1]==3)h=h+6;
		if(a[i-1]==5)h=h+8;
		h++;
	}
	sfile <<"CELLS "<<se<<" "<<h<<endl;
	for (int i=1;i<=se;i++){
			int *element = new int[8];
			a[i-1] = Ng_GetSurfaceElement(mesh,i,element);
			if(a[i-1]==1)sfile <<setprecision(20)<<"3 "<<element[0]<<" "<<element[1]<<" "<<element[2]<<endl;
			if(a[i-1]==2)sfile <<setprecision(20)<<"4 "<<element[0]<<" "<<element[1]<<" "<<element[2]<<" "<<element[3]<<endl;
			if(a[i-1]==3){
				sfile<<"6 ";
				for(int l=0;l<6;l++)sfile<<element[l]<<" ";
				sfile<<endl;
			}
			if(a[i-1]==5){
				sfile<<"8 ";
				for(int l=0;l<8;l++)sfile<<element[l]<<" ";
				sfile<<endl;
			}

			}
	sfile <<endl;
	sfile <<"CELL_TYPES "<< se<<endl;
	for (int i = 0; i < se; i++) {
		if(a[i]==1)sfile<<"5"<<endl;
		if(a[i]==2)sfile<<"9"<<endl;
		if(a[i]==3)sfile<<"22"<<endl;
		if(a[i]==5)sfile<<"23"<<endl;
		}
	sfile.close();

	/*Ng_STL_Uniform_Refinement(stl_geom, mesh);

	std::cout << "elements after refinement: " << Ng_GetNE(mesh) << "\n";
	std::cout << "points   after refinement: " << Ng_GetNP(mesh) << "\n";*/

	Ng_SaveMesh(mesh, "test_ref.vol");
	Ng_Exit();
	return 0;
}












