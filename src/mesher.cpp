
#include "gmsh.h"
#include "tetgen.h"

#include <cstdio>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>

#define TETLIBRARY

using namespace std;


bool meshGMSH(string soubor);
void writeGMSH(string out);

void writeTET(string soubor,tetgenio mesh);

int main(int argc, char **argv)
{
	int opt = 0;
	string infile,outfile,informat,outformat,lib;
	int library=0;
	static struct option long_options[]={
			{"input", required_argument , 0,'i'},
			{"output", required_argument , 0,'o'},
			{"informat", required_argument , 0,'a'},
			{"outformat", required_argument , 0,'b'},
			{"library", required_argument , 0,'l'}
	};
	int long_index=0;
	while((opt = getopt_long(argc, argv, "i:o:a:b:l:",long_options, &long_index))!= -1){
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
		case 'l':{
			cout<< "Library: "<<optarg<<endl;
			lib=optarg;
			break;
				}
		case '?':{
			cout<< "extra: "<<optopt<<endl;
			break;
				}
		}
	}
	infile=infile+"."+informat;
    bool end = true;
    if(lib=="gmsh" || lib=="GMSH")library=1;
    if(lib=="tetgen" || lib=="TETGEN")library=2;
    tetgenio mesh, geo;
    switch (library){
        case 1: {
        	gmsh::initialize();
        	gmsh::option::setNumber("General.Terminal", 1);
        	gmsh::model::add("model");
            end=meshGMSH(infile);
            break;
        }
        case 2: {
        	char file[infile.length()];
        	strcpy(file,infile.c_str());

        	geo.load_stl(file);

       	  	tetgenbehavior tetgen;
       	  	tetgen.object=tetgenbehavior::STL;

       	  	// Settings of mesh
       	  	//tetgen.plc=1;
       	  	//tetgen.quality=1;
       	  	//tetgen.minratio=1.2;
       	  	//tetgen.mindihedral=0;
       	  	tetrahedralize(&tetgen, &geo, &mesh);
       	  	end=true;
       	  	break;
        }
        default:
            cout<<"Wrong chose";
    }

    if (end == true){
    	    string save=outfile+"."+outformat;
    	    switch (library){
            	case 1: {
            		writeGMSH(save);
            		gmsh::finalize();
            		cout<<"Finish"<<endl;
            		break;
            	}
            	case 2: {
            		writeTET(save,mesh);
            		cout<<"Finish"<<endl;
            		break;
            	}
            	default:
            	    cout<<"Wrong chose";
    	    }
    }
	return 0;
}

bool meshGMSH(string soubor){

	  try {
	    gmsh::merge(soubor);
	  } catch(...) {
	    gmsh::logger::write("Could not load file");
	    gmsh::finalize();
	    return false;
	  }
	  double angle = 40;
	  bool forceParametrizablePatches = false;
	  bool includeBoundary = true;
	  double curveAngle = 180;
	  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
	                                      forceParametrizablePatches,
	                                      curveAngle * M_PI / 180.);
	  gmsh::model::mesh::createGeometry();
	  vector<pair<int, int> > s;
	  gmsh::model::getEntities(s, 2);
	  vector<int> sl;
	  for(auto surf : s) sl.push_back(surf.second);
	  int l = gmsh::model::geo::addSurfaceLoop(sl);
	  gmsh::model::geo::addVolume({l});
	  gmsh::model::geo::synchronize();
	  bool funny = true; // false;
	  int f = gmsh::model::mesh::field::add("MathEval");
	  if(funny)
	    gmsh::model::mesh::field::setString(f, "F", "2*Sin((x+y)/5) + 3");
	  else
	    gmsh::model::mesh::field::setString(f, "F", "4");
	  gmsh::model::mesh::field::setAsBackgroundMesh(f);

	  // settings of mesh
	  //gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 10);
	  //gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 0.5);
	  //gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 10);
	  //gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 1.4);
	  //gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.4);
	  //gmsh::model::mesh::setTransfiniteAutomatic();

	  gmsh::model::mesh::generate(3);

	  // Refined mesh
	  //gmsh::model::mesh::recombine();
	  //gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
	  //gmsh::model::mesh::refine();
	  return true;
}

void writeGMSH(string out){
	ofstream sfile;
	sfile.open(out);
	int points=0,element=0,element1=0;
	vector<double> nodecoord;
	vector<pair<int, int>> form;
	vector<size_t> help;

	vector<pair<int, int> > entities;
	gmsh::model::getEntities(entities);
	for(unsigned int i = 0; i < entities.size(); i++) {
		vector<size_t> nodeTags;
	    vector<double> nodeCoords, nodeParams;
	    int dim = entities[i].first, tag = entities[i].second;
	    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

	    points=points+nodeTags.size();
	    for (int k=0;k<nodeCoords.size();k++){
	    	nodecoord.push_back(nodeCoords.at(k));
	    }
	    vector<int> elemTypes;
	    vector<vector<size_t> > elemTags, elemNodeTags;
	    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
	    int numElem = 0;
	    for(unsigned int j = 0; j < elemNodeTags.size(); j++){
	    	for(unsigned int k = 0; k < elemNodeTags[j].size(); k++){
	    		help.push_back(elemNodeTags[j][k]-1);
	    		element1++;
	    	}
	    }
	    for(unsigned int j = 0; j < elemTags.size(); j++){
	          numElem += elemTags[j].size();
	    }
	    element=element+numElem;
	    string type;
	    gmsh::model::getType(dim, tag, type);
	    vector<int> partitions;
	    gmsh::model::getPartitions(dim, tag, partitions);
	    for(unsigned int j = 0; j < elemTypes.size(); j++) {
	    	string name;
	        int d, order, numv, numpv;
	        vector<double> param;
	        gmsh::model::mesh::getElementProperties(elemTypes[j], name, d, order,numv, param, numpv);
	        form.push_back(make_pair(numElem,elemTypes[j]));
	  }
	}
	sfile << "# vtk DataFile Version 2.0"<< endl;
	sfile <<"gmsh, Created by Gmsh"<<endl;
	sfile <<"ASCII"<<endl;
	sfile <<"DATASET UNSTRUCTURED_GRID"<<endl;
	sfile <<"POINTS "<<points<<" double"<<endl;
	for (int k=0;k<nodecoord.size();k++){
		    sfile <<setprecision(20)<< nodecoord.at(k)<< " "<< nodecoord.at(++k)<< " "<< nodecoord.at(++k)<<endl;
	}
	sfile <<endl;
	sfile <<"CELLS "<<element<<" "<<element+element1<<endl;
	int k=0;
	for(unsigned int i=0;i<form.size();i++){
			for(unsigned int j=0; j<form[i].first;j++){

				if(form[i].second==1)sfile<<"2 "<<help[k++]<<" "<<help[k++]<<endl;
				if(form[i].second==2)sfile<<"3 "<<help[k++]<<" "<<help[k++]<<" "<<help[k++]<<endl;
				if(form[i].second==3){
					sfile<<"4 ";
					for(int l=0;l<4;l++)sfile<<help[k++]<<" ";
					sfile<<endl;
				}
				if(form[i].second==4){
					sfile<<"4 ";
					for(int l=0;l<4;l++)sfile<<help[k++]<<" ";
					sfile<<endl;
				}
				if(form[i].second==5){
					sfile<<"8 ";
					for(int l=0;l<8;l++)sfile<<help[k++]<<" ";
					sfile<<endl;
				}
				if(form[i].second==6){
					sfile<<"6 ";
					for(int l=0;l<6;l++)sfile<<help[k++]<<" ";
					sfile<<endl;
				}
				if(form[i].second==7){
					sfile<<"5 ";
					for(int l=0;l<5;l++)sfile<<help[k++]<<" ";
					sfile<<endl;
				}
				if(form[i].second==8)sfile<<"3 "<<help[k++]<<" "<<help[k++]<<" "<<help[k++]<<endl;
				if(form[i].second==9){
									sfile<<"6 ";
									for(int l=0;l<6;l++)sfile<<help[k++]<<" ";
									sfile<<endl;
								}
				if(form[i].second==15)sfile<<"1 "<<help[k++]<<endl;
				if(form[i].second==16){
									sfile<<"8 ";
									for(int l=0;l<8;l++)sfile<<help[k++]<<" ";
									sfile<<endl;
								}
				if(form[i].second==17){
									sfile<<"20 ";
									for(int l=0;l<20;l++)sfile<<help[k++]<<" ";
									sfile<<endl;
								}
			}
		}
	sfile <<endl;
	sfile <<"CELL_TYPES "<< element<<endl;
	for(unsigned int i=0;i<form.size();i++){
		for(unsigned int j=0; j<form[i].first;j++){
			if(form[i].second==1)sfile<<"3"<<endl;
			if(form[i].second==2)sfile<<"5"<<endl;
			if(form[i].second==3)sfile<<"9"<<endl;
			if(form[i].second==4)sfile<<"10"<<endl;
			if(form[i].second==5)sfile<<"12"<<endl;
			if(form[i].second==6)sfile<<"13"<<endl;
			if(form[i].second==7)sfile<<"14"<<endl;
			if(form[i].second==8)sfile<<"21"<<endl;
			if(form[i].second==9)sfile<<"22"<<endl;
			if(form[i].second==15)sfile<<"1"<<endl;
			if(form[i].second==16)sfile<<"23"<<endl;
			if(form[i].second==17)sfile<<"25"<<endl;
		}
	}
	sfile.close();

}

void writeTET(string soubor,tetgenio mesh){
	ofstream sfile;
	sfile.open(soubor);
	sfile << "# vtk DataFile Version 2.0"<< endl;
	sfile <<"tetgen, Created by TETGEN"<<endl;
	sfile <<"ASCII"<<endl;
	sfile <<"DATASET UNSTRUCTURED_GRID"<<endl;
	sfile <<"POINTS " << mesh.numberofpoints << " double"<<endl;
	for (int p = 0; p < mesh.numberofpoints; p++) {
			sfile << mesh.pointlist[p * 3] << " " << mesh.pointlist[p * 3 + 1] << " " << mesh.pointlist[p * 3 + 2] <<endl;
		}
	sfile <<endl;

	sfile << "CELLS " << mesh.numberoftetrahedra << " " << mesh.numberoftetrahedra + 4 * mesh.numberoftetrahedra <<endl;
	for (int e = 0; e < mesh.numberoftetrahedra; e++) {
		sfile << "4";
		for (int en = 0; en < 4; ++en) {
			sfile << " " << mesh.tetrahedronlist[4 * e + en] - 1;
		}
		sfile <<endl;
	}
	sfile <<endl;

	sfile << "CELL_TYPES " << mesh.numberoftetrahedra <<endl;
	for (int e = 0; e < mesh.numberoftetrahedra; e++) {
		sfile << "10"<<endl;
	}
	sfile.close();
}







