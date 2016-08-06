#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>

using namespace plb;
using namespace std;


typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q27Descriptor

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter){
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    
    // Write velocity-norm at x=0.
    Box3D slice(0, nx-1, 0, ny-1, 50, 50);
    imageWriter.writeScaledGif( createFileName("rho", iter, 6),
    *computeDensity(lattice, slice), imSize, imSize);
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter){
        VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
        vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
        vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
}

int main(int argc, char **argv){
    plbInit(&argc, &argv);
    std::string fNameOut = "tmp";

    const plint nx = 100;
    const plint ny = 100;
    const plint nz = 100;
    const T lattice_speed_sound = 1/sqrt(3);

    const T omega = 1.0;
    const T rho0 = 1;
    const T drho = rho0/10;

    const plint maxT = 20000;

    const T lx =  2.3;
    const T dx =lx/nx;

    Array<T,3> u0(0, 0, 0);

    plint cx = util::roundToInt(nx/2);
    plint cy = util::roundToInt(ny/2);
    plint cz = util::roundToInt(nz/2);
    Array<T,3> centerLB(cx , cy, cz);

    global::directories().setOutputDir(fNameOut+"/");

    pcout << "Creation of the lattice." << endl;
    MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx,ny,nz, new CompleteBGKdynamics<T,DESCRIPTOR>(omega));

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), rho0 , u0 );


    AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
    new AnechoicDynamics<T,DESCRIPTOR>(omega);
    T delta_efective = 30;
    anechoicDynamics->setDelta(delta_efective);
    anechoicDynamics->setRhoBar_target(0);
    anechoicDynamics->setJ_target(Array<T,3>(0, 0, 0));
    //DotList3D points_to_aplly_dynamics;
    //points_to_aplly_dynamics.addDot(Dot3D(100,100,100));
    defineDynamics(lattice, Box3D( 10, 30 , 50, 70 , 40,70), anechoicDynamics);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;
    pcout << std::endl << "dx:" << dx << std::endl;

    pcout << "Simulation begins" << endl;
    for (plint iT=0; iT<maxT; ++iT){
            
        if (iT== 0){
            Box3D impulse( 50, 50 , 50, 50 , 50,50);
            initializeAtEquilibrium( lattice, impulse, rho0+drho*5 , u0 );
            lattice.initialize();

        }

        if (iT % 20 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
            writeGifs(lattice,iT);
            writeVTK(lattice, iT);
        }

        lattice.collideAndStream();

    }

    pcout << "End of simulation at iteration " << endl;

}
