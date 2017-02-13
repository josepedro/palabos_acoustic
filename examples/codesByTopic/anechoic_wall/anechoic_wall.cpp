#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <time.h>
#include <cfloat> 

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR MRTD3Q19Descriptor
typedef MRTdynamics<T,DESCRIPTOR> BackgroundDynamics;
typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;

// ---------------------------------------------
// Includes of acoustics resources
#include "acoustics/acoustics3D.h"
using namespace plb_acoustics_3D;
// ---------------------------------------------

int main(int argc, char **argv){
    plbInit(&argc, &argv);

    global::timer("mainLoop").start();

    const T rho0 = 1;
    const T drho = rho0/100;
    const plint nx = 100;
    const plint ny = 100;
    const plint nz = 250;
    const T omega = 1.985;
    const plint maxT = 10000;
    const T lattice_speed_sound = 1/sqrt(3);
    const T cs2 = lattice_speed_sound*lattice_speed_sound;
    Array<T,3> u0(0, 0, 0);

    // Saving a propery directory
    std::string fNameOut = currentDateTime() + "+tmp";
    std::string command = "mkdir -p " + fNameOut;
    char to_char_command[1024];
    strcpy(to_char_command, command.c_str());
    system(to_char_command);
    global::directories().setOutputDir(fNameOut+"/");


    // Setting anechoic dynamics like this way
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(nx, ny, nz,  new AnechoicBackgroundDynamics(omega));
    defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));

    pcout << "Creation of the lattice. Time max: " << maxT << endl;

    pcout << getMultiBlockInfo(lattice) << std::endl;

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0 , u0);

    T rhoBar_target = 0;
    Array<T,3> j_target(0, 0, 0);
    T size_anechoic_buffer = 30;
    defineAnechoicMRTBoards(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;

    pcout << "Simulation begins" << endl;

    // Setting probes to evaluate coefficient reflection of ABC
    Array<plint,3> position_A(nx/2, ny/2, nz - size_anechoic_buffer);
    Array<plint,3> position_B(nx/2, ny - size_anechoic_buffer, nz - size_anechoic_buffer);
    Array<plint,3> position_C(nx/2, ny - size_anechoic_buffer, nz/2);
    Array<plint,3> position_D(nx/2, 0.75*ny, size_anechoic_buffer + 1);
    string name_abc = "abc";
    Coefficient_Reflection_Probes coefficient_reflection_probes(position_A,
    position_B, position_C, position_D, fNameOut, name_abc);
    // ---------------------------------------------------------

    // Recording entering signal -------------------------------
    std::string signal_in_string = fNameOut+"/signal_in.dat";
    char to_char_signal_in[1024];
    strcpy(to_char_signal_in, signal_in_string.c_str());
    plb_ofstream history_signal_in(to_char_signal_in);
    // ---------------------------------------------------------

    // Mean for-loop
    for (plint iT=0; iT<maxT; ++iT){
        if (iT <= maxT){
            T freq_omega_max = 2*M_PI*0.4;
            plint total_signals = 20;
	        T chirp_hand = get_linear_chirp_AZ_freq_omega(freq_omega_max, 
                total_signals, maxT, iT, drho);
            //T chirp_hand = get_linear_chirp(ka_min, ka_max, maxT_final_source, iT, drho, radius);
            //T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            history_signal_in << setprecision(10) << chirp_hand << endl;
            Array<T,3> j_target(0, 0, 0);
            Box3D signal_in(nx/2, nx/2, ny/2, ny/2, nz/4, nz/4);
            initializeAtEquilibrium(lattice, signal_in, chirp_hand, j_target);
        }

        if (iT % 50 == 0) {
            pcout << "Iteration " << iT << endl;
             pcout << " energy ="
            << setprecision(10) << getStoredAverageEnergy<T>(lattice)
            << " rho ="
            << getStoredAverageDensity<T>(lattice)
            << " max_velocity ="
            << setprecision(10) << (getStoredMaxVelocity<T>(lattice))/lattice_speed_sound
            << endl;
        }

        if (iT == (maxT/2) + (maxT/4)) {
            //writeGifs(lattice,iT);
            //writeVTK(lattice, iT, rho0, drho);
        }

        if (iT%50 == 0) {
            //writeGifs(lattice,iT);
            //writeVTK(lattice, iT, rho0, drho);
        }


        coefficient_reflection_probes.save_point(lattice, rho0, cs2);
        lattice.collideAndStream();
    }

    T total_time_simulation = global::timer("mainLoop").stop();
    pcout << "End of simulation at iteration with total time: " << total_time_simulation << endl;
}
