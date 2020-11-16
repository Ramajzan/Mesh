
#include "gmsh.h"

#include <cstdio>
#include <iostream>
using namespace std;

bool meshGMSH(string soubor);


int main(int argc, char **argv)
{
	// parse args (e.g. https://www.gnu.org/software/libc/manual/html_node/Getopt.html)
	cout << "Nacti soubor" << std::endl;
	string soubor;
	cin>>soubor;
	cout<< "Soubor: "<< soubor << std::endl;
    cout<< "Jakou chceš použít knihovnu: 1) GMSH ";
    int a;
    bool end = true;
    cin>> a;
    switch (a){
        case 1: {
        	gmsh::initialize();
        	gmsh::option::setNumber("General.Terminal", 1);
        	gmsh::model::add("t13");
            end=meshGMSH(soubor);
            break;
        }
        default:
            cout<<"Chybný výběr";
    }

    if (end == true){
    		string format;
    		cout<<"Vyber formát: "<<endl;
    	    cin>>format;
    	    int help=soubor.find(".");
    	    soubor=soubor.substr(0,help+1);
    	    soubor=soubor+format;
    	    switch (a){
            	case 1: {
            		gmsh::write(soubor);
            		gmsh::finalize();
            		break;
            	}
    	    }
    }

	// store to selected format (manually)

	return 0;
}

bool meshGMSH(std::string soubor){

	  try {
	    gmsh::merge(soubor);
	  } catch(...) {
	    gmsh::logger::write("Could not load STL mesh: bye!");
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
	  std::vector<std::pair<int, int> > s;
	  gmsh::model::getEntities(s, 2);
	  std::vector<int> sl;
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
	  gmsh::model::mesh::generate(3);

	  return true;
}









