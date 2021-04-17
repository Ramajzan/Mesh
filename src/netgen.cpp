
#include <cstddef>
#define OCCGEOMETRY

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

#include <ctime>

using namespace std;
using namespace nglib;

static int check(Ng_Result result)
{
	switch (result) {
	case NG_ERROR: cout << "NG ERROR\n"; break;
	case NG_OK: break;
	case NG_SURFACE_INPUT_ERROR: cout << "NG SURFACE INPUT ERROR\n"; break;
	case NG_VOLUME_FAILURE: cout << "NG VOLUME FAILURE\n"; break;
	case NG_STL_INPUT_ERROR: cout << "NG STL INPUT ERROR\n"; break;
	case NG_SURFACE_FAILURE: cout << "NG SURFACE FAILURE\n"; break;
	case NG_FILE_NOT_FOUND: cout << "NG FILE NOT FOUND\n"; break;
	}
	return result != NG_OK;
}

int main(int argc, char **argv)
{
	clock_t start;
	double duration;
	int opt = 0;
	string infile,outfile,informat,outformat,lib;
	int library=0;
	static struct option long_options[]={
			{"input", required_argument , 0,'i'},
			{"output", required_argument , 0,'o'},
			{"informat", required_argument , 0,'a'},
			{"outformat", required_argument , 0,'b'}
	};
	int long_index=0;
	while((opt = getopt_long(argc, argv, "i:o:a:b:",long_options, &long_index))!= -1){
		switch(opt){
			case 'i':{
				cout<< "Input file: "<<optarg<<endl;
				infile=optarg;
				break;
				}
			case 'o':{
				cout<< "Output file: "<<optarg<<endl;
				outfile=optarg;
				break;
				}
			case 'a':{
				cout<< "Input format: "<<optarg<<endl;
				informat=optarg;
				break;
				}
			case 'b':{
				cout<< "Output format: "<<optarg<<endl;
				outformat=optarg;
				break;
				}
			case '?':{
				cout<< "extra: "<<optopt<<endl;
				break;
				}
		}
	}
	infile=infile+"."+informat;
	char file[infile.length()+1];
	strcpy(file,infile.c_str());

	Ng_Meshing_Parameters mp;

		//!< Switch to enable / disable usage of local mesh size modifiers
	mp.uselocalh=true;
		//!< Maximum global mesh size allowed
	mp.maxh;
		//!< Minimum global mesh size allowed
	mp.minh;
		//!< Mesh density: 0...1 (0 => coarse; 1 => fine)
	mp.fineness;
		//!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
	mp.grading;

		//!< Number of elements to generate per edge of the geometry
	mp.elementsperedge;
		//!< Elements to generate per curvature radius
	mp.elementspercurve;

		//!< Enable / Disable mesh refinement at close edges
	mp.closeedgeenable;
		//!< Factor to use for refinement at close edges (larger => finer)
	mp.closeedgefact;

		//!< Enable / Disable user defined minimum edge length for edge subdivision
	mp.minedgelenenable;
		//!< Minimum edge length to use while subdividing the edges (default = 1e-4)
	mp.minedgelen;

		//!< Generate second-order surface and volume elements
	mp.second_order;
		//!< Creates a Quad-dominated mesh
	mp.quad_dominated;

		//!< Optional external mesh size file
	mp.meshsize_filename;

		//!< Enable / Disable automatic surface mesh optimization
	mp.optsurfmeshenable;
		//!< Enable / Disable automatic volume mesh optimization
	mp.optvolmeshenable;

		//!< Number of optimize steps to use for 3-D mesh optimization
	mp.optsteps_3d;
		//!< Number of optimize steps to use for 2-D mesh optimization
	mp.optsteps_2d;

		// Philippose - 13/09/2010
		// Added a couple more parameters into the meshing parameters list
		// from Netgen into Nglib
		//!< Invert all the volume elements
	mp.invert_tets;
		//!< Invert all the surface triangle elements
	mp.invert_trigs;

		//!< Check for overlapping surfaces during Surface meshing
	mp.check_overlap;
		//!< Check for overlapping surface elements before volume meshing
	mp.check_overlapping_boundary;

    start=clock();
	Ng_Init();
	Ng_Mesh *mesh = Ng_NewMesh();
	mp.Transfer_Parameters();

	if(informat=="stl"){
		Ng_STL_Geometry *geom = Ng_STL_LoadGeometry(file);
		check(Ng_STL_InitSTLGeometry(geom));
		check(Ng_STL_MakeEdges(geom, mesh, &mp));
		check(Ng_STL_GenerateSurfaceMesh(geom, mesh, &mp));
		check(Ng_GenerateVolumeMesh(mesh, &mp));
	}
	if(informat=="stp"){
		Ng_OCC_Geometry *geom = Ng_OCC_Load_STEP(file);
		check(Ng_OCC_SetLocalMeshSize(geom, mesh, &mp));
		check(Ng_OCC_GenerateEdgeMesh(geom, mesh, &mp));
		check(Ng_OCC_GenerateSurfaceMesh(geom, mesh, &mp));
		check(Ng_GenerateVolumeMesh(mesh, &mp));
	}
	if(informat=="igs"){
		Ng_OCC_Geometry *geom = Ng_OCC_Load_IGES(file);
		check(Ng_OCC_SetLocalMeshSize(geom, mesh, &mp));
		check(Ng_OCC_GenerateEdgeMesh(geom, mesh, &mp));
		check(Ng_OCC_GenerateSurfaceMesh(geom, mesh, &mp));
		check(Ng_GenerateVolumeMesh(mesh, &mp));
		}
	int points = Ng_GetNP(mesh), elements = Ng_GetNE(mesh);

	string save=outfile+"."+outformat;
	ofstream sfile;
	sfile.open(save);
	sfile << "# vtk DataFile Version 2.0" << endl;
	sfile <<"NGsolve, Created by NGsolve" << endl;
	sfile <<"ASCII" << endl;
	sfile <<"DATASET UNSTRUCTURED_GRID" << endl;
	sfile <<"POINTS " << points << " double" << endl;

	double point[3];
	for (int p = 1; p <= points; p++) {
			Ng_GetPoint(mesh, p, point);
			sfile << point[0] << " " << point[1] << " " << point[2] << endl;
		}
	sfile << endl;

	int element[NG_VOLUME_ELEMENT_MAXPOINTS], enodes = 0;
	for (int e = 1; e <= elements; e++) {
		switch (Ng_GetVolumeElement (mesh, e, element)) {
			case NG_TET: enodes += 4; break;
			case NG_PYRAMID: enodes += 5; break;
			case NG_PRISM: enodes += 6; break;
			case NG_TET10: enodes += 10; break;
			}
	}

	sfile << "CELLS " << elements << " " << elements + enodes << endl;
	for (int e = 1; e <= elements; e++) {
		switch (Ng_GetVolumeElement (mesh, e, element)) {
			case NG_TET: enodes = 4; break;
			case NG_PYRAMID: enodes = 5; break;
			case NG_PRISM: enodes = 6; break;
			case NG_TET10: enodes = 10; break;
		}
		sfile << enodes;
		for (int en = 0; en < enodes; ++en) {
			sfile << " " << element[en] - 1;
		}
		sfile << endl;
	}
	sfile << endl;

	sfile << "CELL_TYPES " << elements << endl;
	for (int e = 1; e <= elements; e++) {
		switch (Ng_GetVolumeElement (mesh, e, element)) {
			case NG_TET: sfile << "10" << endl; break;
			case NG_PYRAMID: sfile << "14" << endl;; break;
			case NG_PRISM: sfile << "13" << endl; break;
			case NG_TET10: sfile << "24" << endl; break;
		}
	}
	cout<<"nodes: "<<points<<"\n elements: "<<elements<<endl;
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Doba meshe: "<<duration<<endl;
	return 0;

}












