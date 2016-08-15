#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>

using namespace plb;
using namespace std;


typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q27Descriptor

// ---------------------------------------------
// Includes of acoustics resources
#include "acoustics/acoustics3D.h"
using namespace plb_acoustics_3D;
// ---------------------------------------------

const T rho0 = 1;
const T drho = rho0/10;

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter){
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    
    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    imageWriter.writeScaledGif( createFileName("velocity", iter, 6),
    *computeDensity(lattice, slice), imSize, imSize);

    imageWriter.writeGif( createFileName("rho", iter, 6),
    *computeDensity(lattice, slice), 
    (T) rho0 - drho/100000, (T) rho0 + drho/100000, imSize, imSize);
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

    const T omega = 1.9;
    const plint maxT = 120000;

    const T lx =  2.3;
    const T dx =lx/nx;

    Array<T,3> u0(0, 0, 0);

    plint cx = util::roundToInt(nx/2);
    plint cy = util::roundToInt(ny/2);
    plint cz = util::roundToInt(nz/2);
    Array<T,3> centerLB(cx , cy, cz);

    global::directories().setOutputDir(fNameOut+"/");

    pcout << "Creation of the lattice." << endl;
    MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx,ny,nz, 
        new CompleteBGKdynamics<T,DESCRIPTOR>(omega));

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), rho0 , u0 );

    plint size_square = 4;
    Box3D square(
    nx/2 - size_square/2, nx/2 + size_square/2,
    ny/2 - size_square/2, ny/2 + size_square/2, 
    nz/2 - size_square/2, nz/2 + size_square/2);
    defineDynamics(lattice, square, new BounceBack<T,DESCRIPTOR>((T)0));

    // Anechoic Condition
    T rhoBar_target = 0;
    const T mach_number = 0.2;
    const T velocity_flow = mach_number*lattice_speed_sound;
    Array<T,3> j_target(velocity_flow, 0.0/std::sqrt(3), 0.0/std::sqrt(3));
    T size_anechoic_buffer = 30;
    // Define Anechoic Boards
    defineAnechoicBoards(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target);
    

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;
    pcout << std::endl << "dx:" << dx << std::endl;

    pcout << "Simulation begins" << endl;

    plb_ofstream history_pressures("history_pressures.dat");
    for (plint iT=0; iT<maxT; ++iT){
        if (iT == 0){
            //Box3D impulse(nx/2, nx/2, ny/2, ny/2, nz/2, nz/2);
            //initializeAtEquilibrium( lattice, impulse, rho0 + drho, u0 );
        }

        if (iT % 100 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
              pcout << " energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << " rho ="
                  << getStoredAverageDensity<T>(lattice)
                  << " max_velocity ="
                  << setprecision(10) << (getStoredMaxVelocity<T>(lattice))/lattice_speed_sound
                  << endl;
            writeGifs(lattice,iT);
            //writeVTK(lattice, iT);
        }


        history_pressures << setprecision(10) << lattice.get(nx/2+30, ny/2+30, nz/2+30).computeDensity() - rho0 << endl;
        lattice.collideAndStream();

    }

    pcout << "End of simulation at iteration " << endl;

}
