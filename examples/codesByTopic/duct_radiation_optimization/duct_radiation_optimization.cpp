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

    const T rho0 = 1;
    const T drho = rho0/100;
    const plint radius = 20;
    const plint diameter = 2*radius;
    const plint nx = 6*diameter + 60;
    //const plint nx = 2*diameter + 60;
    const plint ny = 6*diameter + 60;
    //const plint ny = 2*diameter + 60;
    const plint position_duct_z = 0;
    const plint length_duct = 30 + 120 + 5*113 + 3*diameter;
    //const plint length_duct = 3*diameter;
    const plint nz = length_duct + 2*3*diameter + 30;
    //const plint nz = length_duct + 3*diameter + 30;
    const T lattice_speed_sound = 1/sqrt(3);
    const T omega = 1.985;
    const plint maxT = 2*(pow(2,13) + nz*sqrt(3));
    Array<T,3> u0(0, 0, 0);
    const Array<plint,3> position(nx/2, ny/2, position_duct_z);
    const plint thickness_duct = 2;
    const plint radius_intern = radius - 2;
    const plint maxT_final_source = maxT - nz*sqrt(3);
    const T ka_max = 2.5;
    const T ka_min = 0;
    const T cs2 = lattice_speed_sound*lattice_speed_sound;
    clock_t t;

    // Saving a propery directory
    std::string fNameOut = currentDateTime() + "+tmp";
    std::string command = "mkdir -p " + fNameOut;
    char to_char_command[1024];
    strcpy(to_char_command, command.c_str());
    system(to_char_command);
    std::string command_copy_script_matlab = "cp duct_radiation.m " + fNameOut;
    char to_char_command_copy_script_matlab[1024];
    strcpy(to_char_command_copy_script_matlab, command_copy_script_matlab.c_str());
    system(to_char_command_copy_script_matlab);
    global::directories().setOutputDir(fNameOut+"/");


    // Setting anechoic dynamics like this way
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(nx, ny, nz,  new AnechoicBackgroundDynamics(omega));
    defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));

    pcout << "Creation of the lattice. Time max: " << maxT << endl;
    pcout << "Duct radius: " << radius << endl;

    pcout << getMultiBlockInfo(lattice) << std::endl;

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0 , u0);

    // Set NoDynamics to improve performance!
    plint off_set_z = position_duct_z + length_duct - 3*diameter - 30;
    set_nodynamics(lattice, nx, ny, off_set_z);
        
    T rhoBar_target = 0;
    const T mach_number = 0.2;
    const T velocity_flow = mach_number*lattice_speed_sound;
    Array<T,3> j_target(0, 0, 0);
    T size_anechoic_buffer = 30;
    defineAnechoicMRTBoards_limited(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target, off_set_z);

    //const T mach_number = 0;
    build_duct(lattice, nx, ny, position, radius, length_duct, thickness_duct, omega);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;

    pcout << "Simulation begins" << endl;

    // Setting probes ------------------------------------------
    plint distance_group_A = 28;
    plint distance_group_B = 113 + 123;
    System_Abom_Measurement system_abom_measurement(lattice, position, radius,
    distance_group_A, distance_group_B, fNameOut);

    plint double_microphone_distance = 5;

    plint distance_boca = 0*radius;
    string name_boca = "boca";
    Two_Microphones two_microphones_boca(radius, double_microphone_distance, 
            length_duct, position, fNameOut, name_boca, distance_boca, nx, ny, nz);

    plint distance_1r = 1*radius;
    string name_1r = "1r";
    Two_Microphones two_microphones_1r(radius, double_microphone_distance, 
            length_duct, position, fNameOut, name_1r, distance_1r, nx, ny, nz);

    plint distance_2r = 2*radius;
    string name_2r = "2r";
    Two_Microphones two_microphones_2r(radius, double_microphone_distance, 
            length_duct, position, fNameOut, name_2r, distance_2r, nx, ny, nz);

    plint distance_3r = 3*radius;
    string name_3r = "3r";
    Two_Microphones two_microphones_3r(radius, double_microphone_distance, 
            length_duct, position, fNameOut, name_3r, distance_3r, nx, ny, nz);

    plint distance_6r = 6*radius;
    string name_6r = "6r";
    Two_Microphones two_microphones_6r(radius, double_microphone_distance, 
            length_duct, position, fNameOut, name_6r, distance_6r, nx, ny, nz);

    // ---------------------------------------------------------
    

    // Recording entering signal -------------------------------
    std::string signal_in_string = fNameOut+"/signal_in.dat";
    char to_char_signal_in[1024];
    strcpy(to_char_signal_in, signal_in_string.c_str());
    plb_ofstream history_signal_in(to_char_signal_in);
    // ---------------------------------------------------------

    // Important information about simulation ------------------    
    t = clock();
    std::string AllSimulationInfo_string = fNameOut + "/AllSimulationInfo.txt";
    char to_char_AllSimulationInfo[1024];
    strcpy(to_char_AllSimulationInfo, AllSimulationInfo_string.c_str());
    plb_ofstream AllSimulationInfo(to_char_AllSimulationInfo);
    
    std::string title = "\nRODANDO AGORA COM ESCOAMENTO E A CAVIDADE MAIOR (2 VESZES PARA FRENTE).\n"; 
    
    AllSimulationInfo << endl
    << title << endl
    << "Dados da simulação" << endl
    << "Lattice:" << endl << endl 
    << "nx: " << nx << " ny: " << ny << " nz: " << nz << endl
    << " omega: " << omega << endl << endl
    << "Tempos: " << endl
    << "Total Time step: " << maxT << endl
    << "Raio do duto: " << radius << endl
    << "Espessura: " << thickness_duct << endl
    << "Tamanho duto: " << length_duct << endl
    << "Posicao do duto: " << position[2] << endl;
    // --------------------------------------------------------


    //pcout << "!!Loading lattice initial condition!!" << endl;
    //loadBinaryBlock(lattice, "checkpoint.dat");
    // Mean for-loop
    for (plint iT=0; iT<maxT; ++iT){
        if (iT <= maxT_final_source && iT > maxT/2){
            plint total_signals = 20;
	    T chirp_hand = get_linear_chirp_AZ(ka_max,  total_signals, maxT_final_source, iT - maxT/2, drho, radius);
            //T chirp_hand = get_linear_chirp(ka_min, ka_max, maxT_final_source, iT, drho, radius);
            //T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            history_signal_in << setprecision(10) << chirp_hand << endl;
            Array<T,3> j_target(0, 0, velocity_flow);
            set_source(lattice, position, chirp_hand, j_target, radius, radius_intern, nx, ny);
        }else{
              Array<T,3> j_target(0, 0, velocity_flow);
            set_source(lattice, position, rho0, j_target, radius, radius_intern, nx, ny);
        }

        if (iT % 100 == 0) {
            pcout << "Iteration " << iT << endl;
             pcout << " energy ="
            << setprecision(10) << getStoredAverageEnergy<T>(lattice)
            << " rho ="
            << getStoredAverageDensity<T>(lattice)
            << " max_velocity ="
            << setprecision(10) << (getStoredMaxVelocity<T>(lattice))/lattice_speed_sound
            << endl;
        }

        if (iT % 200 == 0) {
            //writeGifs(lattice,iT);
            //writeVTK(lattice, iT, rho0, drho);
        }

        if (iT == maxT/2){
            pcout << "Saving the state of the simulation ..." << endl;
            //saveRawMultiBlock(lattice, "checkpoint.dat");
            saveBinaryBlock(lattice, "checkpoint.dat");
        }

        // extract values of pressure and velocities
        system_abom_measurement.save_point(lattice, rho0, cs2);
        two_microphones_1r.save_point(lattice, rho0, cs2);
        two_microphones_2r.save_point(lattice, rho0, cs2);
        two_microphones_3r.save_point(lattice, rho0, cs2);
        two_microphones_6r.save_point(lattice, rho0, cs2);
        two_microphones_boca.save_point(lattice, rho0, cs2);

        lattice.collideAndStream();
    }

    t = (clock() - t)/CLOCKS_PER_SEC;
    AllSimulationInfo << endl << "Execution time: " << t << " segundos" << endl;

    pcout << "End of simulation at iteration " << endl;
}
