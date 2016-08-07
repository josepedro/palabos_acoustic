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
const T drho = rho0*7;

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

    const plint nx = 100;
    const plint ny = 100;
    const plint nz = 100;
    const T lattice_speed_sound = 1/sqrt(3);

    const T omega = 1.0;
    const plint maxT = 210*2;

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

    // Anechoic Condition
    T rhoBar_target = 0;
    Array<T,3> j_target(0, 0, 0);
    T size_anechoic_buffer = 30;
    // Define Anechoic Boards
    defineAnechoicBoards(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target);
    

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;
    pcout << std::endl << "dx:" << dx << std::endl;

    pcout << "Simulation begins" << endl;
    for (plint iT=0; iT<maxT; ++iT){
        T lattice_speed_sound = 1/sqrt(3);
        T rho_changing = rho0 + drho*sin(2*M_PI*4*(lattice_speed_sound/20)*iT);        
        Box3D impulse(nx/2, nx/2, ny/2, ny/2, nz/2, nz/2);
        initializeAtEquilibrium( lattice, impulse, rho_changing, u0 );


        if (iT % 5 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
            writeGifs(lattice,iT);
            //writeVTK(lattice, iT);
        }

        if (iT == 200){
            plb_ofstream sinusoidal_001_points("tmp/sinusoidal_001_points.dat");
            plb_ofstream sinusoidal_010_points("tmp/sinusoidal_010_points.dat");
            plb_ofstream sinusoidal_011_points("tmp/sinusoidal_011_points.dat");
            plb_ofstream sinusoidal_100_points("tmp/sinusoidal_100_points.dat");
            plb_ofstream sinusoidal_101_points("tmp/sinusoidal_101_points.dat");
            plb_ofstream sinusoidal_110_points("tmp/sinusoidal_110_points.dat");
            plb_ofstream sinusoidal_111_points("tmp/sinusoidal_111_points.dat");

            for (plint z = 0; z < nz; ++z){
                for (plint y = 0; y < ny; ++y){
                    for (plint x = 0; x < nx; ++x){
                        if (x == 0 && y == 0){
                            sinusoidal_001_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (x == 0 && z == 0){
                            sinusoidal_010_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (x == 0 && y == z){
                            sinusoidal_011_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (y == 0 && z == 0){
                            sinusoidal_100_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (x == z && y == 0){
                            sinusoidal_101_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (x == y && z == 0){
                            sinusoidal_110_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                        if (x == y && z == x){
                            sinusoidal_111_points << setprecision(10) <<
                            lattice.get(x,y,z).computeDensity() - rho0 << endl;
                        }
                    }
                }                
            }
        }

        lattice.collideAndStream();

    }

    pcout << "End of simulation at iteration " << endl;

}
