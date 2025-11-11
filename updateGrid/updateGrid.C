#include "fvCFD.H"
#include <array>
#include <vector>
#include <cmath>

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "functionObject.H"

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    std::string inFileName = "cfd_grid.csv";
    int tot = 0;
    std::string line;
    std::ifstream myfile(inFileName);
    while (std::getline(myfile, line))
        ++tot;
    
    std::ifstream inFile;
    inFile.open(inFileName.c_str());
    scalar x, y, z = 0;
    std::ofstream outdata;
    outdata.open("cfdOutput.csv");
    for (int c=0; c<tot; c++) {
	inFile >> x >> y >> z;
 
	vector pos = {x,y,z};
        const label celli = mesh.findCell(pos);

	if (celli == -1) {
            outdata << 0 << " " << 0 << " " << 0 << std::endl;
	} else {
	    vector vel = U[celli];
	    outdata << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;
	}
    }
    outdata.close();

    return 0;
}
