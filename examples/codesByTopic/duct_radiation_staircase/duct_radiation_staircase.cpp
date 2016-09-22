#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>

using namespace plb;
using namespace plb::descriptors;
using namespace std;


typedef double T;
typedef Array<T,3> Velocity;
//#define DESCRIPTOR descriptors::D3Q27Descriptor
#define DESCRIPTOR MRTD3Q19Descriptor
 typedef MRTdynamics<T,DESCRIPTOR> BackgroundDynamics;
 typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;

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
    //imageWriter.writeGif(createFileName("u", iT, 6),
     //*computeDensity(lattice), );

    imageWriter.writeGif( createFileName("rho", iter, 6),
    *computeDensity(lattice, slice), 
    (T) rho0 - drho/1000000, (T) rho0 + drho/1000000, imSize, imSize);
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter){
        VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
        vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
        vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
}

int main(int argc, char **argv){
    plbInit(&argc, &argv);
    std::string fNameOut = "tmp";

    const plint length_domain = 420;
    const plint nx = length_domain;
    const plint ny = length_domain;
    const plint nz = length_domain;
    const T lattice_speed_sound = 1/sqrt(3);

    const T omega = 1.985;
    const plint maxT = 120000;

    Array<T,3> u0(0, 0, 0);

    global::directories().setOutputDir(fNameOut+"/");

    // Setting anechoic dynamics like this way
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(nx, ny, nz,  new AnechoicBackgroundDynamics(omega));
    defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));

    pcout << "Creation of the lattice." << endl;

    pcout << getMultiBlockInfo(lattice) << std::endl;

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), rho0 , u0 );

    plint size_square = 10;
    Box3D square(
    35, 35 + size_square,
    ny/2 - size_square/2, ny/2 + size_square/2, 
    nz/2 - size_square/2, nz/2 + size_square/2);
    defineDynamics(lattice, square, new BounceBack<T,DESCRIPTOR>((T)0));
    
    T rhoBar_target = 0;
    const T mach_number = 0.2;
    const T velocity_flow = mach_number*lattice_speed_sound;
    Array<T,3> j_target(velocity_flow, 0, 0);
    T size_anechoic_buffer = 20;
    defineAnechoicMRTBoards(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;

    pcout << "Simulation begins" << endl;

    plb_ofstream history_pressures("tmp/history_pressures.dat");
    plb_ofstream history_velocities_x("tmp/history_velocities_x.dat");
    plb_ofstream history_velocities_y("tmp/history_velocities_y.dat");
    plb_ofstream history_velocities_z("tmp/history_velocities_z.dat");
    for (plint iT=0; iT<maxT; ++iT){
        if (iT != 0){
            T lattice_speed_sound = 1/sqrt(3);
            T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            Box3D impulse(nx/2 + 20, nx/2 + 20, ny/2 + 20, ny/2 + 20, nz/2 + 20, nz/2 + 20);
            //initializeAtEquilibrium( lattice, impulse, rho_changing, u0 );
        }

        if (iT % 10 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
            //writeGifs(lattice,iT);
            writeVTK(lattice, iT);
        }

        history_pressures << setprecision(10) << lattice.get(nx/2+30, ny/2+30, nz/2+30).computeDensity() - rho0 << endl;
        Array<T,3> velocities;
        lattice.get(nx/2+30, ny/2+30, nz/2+30).computeVelocity(velocities);
        history_velocities_x << setprecision(10) << velocities[0]/lattice_speed_sound << endl;
        history_velocities_y << setprecision(10) << velocities[1]/lattice_speed_sound << endl;
        history_velocities_z << setprecision(10) << velocities[2]/lattice_speed_sound << endl;
        lattice.collideAndStream();

    }

    pcout << "End of simulation at iteration " << endl;
}
