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

const T rho0 = 1;
const T drho = rho0/100;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

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
        std::auto_ptr<MultiScalarField3D<T> > velocity(plb::computeVelocityComponent(lattice, lattice.getBoundingBox(), 2));
        
        vtkOut.writeData<T>(*velocity, "velocity", 1.);
}

void build_duct(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint nx, plint ny,
    Array<plint,3> position, plint radius, plint length, plint thickness, T omega){
    length += 4;
    plint anechoic_size = 20;
    // Duct is constructed along the Z direction
    //plint size_square = 50;
    plint size_square = 2*radius;
    plint radius_intern = radius - thickness;
    for (plint x = position[0] - radius; x < nx/2 + size_square/2; ++x){
        for (plint y = position[1] - radius; y < ny/2 + size_square/2; ++y){
            for (plint z = position[2]; z < length + position[2]; ++z){

                if (radius*radius > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2)){
                    DotList3D points_to_aplly_dynamics;
                    points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                    defineDynamics(lattice, points_to_aplly_dynamics, new BounceBack<T,DESCRIPTOR>(0));
                }
                // extrude
                if (radius_intern*radius_intern > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2) && z > position[2] + 2){
                    DotList3D points_to_aplly_dynamics;
                    points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                    if (z < position[2] + 2 + anechoic_size && z > position[2] + 3){
                        Array<T,3> u0(0, 0, 0);
                        T rhoBar_target = 0;
                        AnechoicBackgroundDynamics *anechoicDynamics = 
                        new AnechoicBackgroundDynamics(omega);
                        T delta_efective = anechoic_size - z - position[2] - 2;
                        anechoicDynamics->setDelta(delta_efective);
                        anechoicDynamics->setRhoBar_target(rhoBar_target);
                        anechoicDynamics->setJ_target(u0);
                        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        
                    }else{
                        defineDynamics(lattice, points_to_aplly_dynamics, new BackgroundDynamics(omega));
                    }
                }
            }
        }
    }
}

T get_linear_chirp_value(T ka_min, T ka_max, plint maxT_final_source, plint iT, T drho, T radius){
    T lattice_speed_sound = 1/sqrt(3);
    T initial_frequency = ka_min*lattice_speed_sound/(2*M_PI*radius);
    T frequency_max_lattice = ka_max*lattice_speed_sound/(2*M_PI*radius);
    T variation_frequency = (frequency_max_lattice - initial_frequency)/maxT_final_source;
    T frequency_function = initial_frequency*iT + (variation_frequency*iT*iT)/2;
    T phase = 2*M_PI*frequency_function;
    T chirp_hand = 1. + drho*sin(phase);

    return chirp_hand;
}

T get_linear_chirp_AZ(T ka_max, plint total_signals, plint maxT_final_source, plint iT, T drho, T radius){
    T cs = 1/sqrt(3);
    T chirp_value = 1;
    for (plint n_signal = 1; n_signal <= total_signals; n_signal++){
        T interval = ka_max/total_signals;
        T phase = (n_signal*interval*cs*iT)/radius;
        chirp_value += drho*sin(phase);
    }
    return chirp_value;
}

int main(int argc, char **argv){
    plbInit(&argc, &argv);

    //const plint length_domain = 420;
    const plint radius = 50;
    const plint diameter = 2*radius;
    //const plint length_domain = 150;
    const plint nx = 6*diameter + 60;
    const plint ny = 6*diameter + 60;
    const plint length_duct = 6*diameter + 20;
    const plint position_duct_z = 30;
    const plint nz = length_duct + 3*diameter + 60;
    const T lattice_speed_sound = 1/sqrt(3);
    const T omega = 1.985;
    const plint maxT = pow(2,13) + nz*sqrt(3);
    Array<T,3> u0(0, 0, 0);

    //const plint maxT = 2*120/lattice_speed_sound;
    const plint maxT_final_source = maxT - nz*sqrt(3);
    const T ka_max = 1.86;
    //const T ka_min = 0;
    const T cs2 = lattice_speed_sound*lattice_speed_sound;
    const plint source_radius = radius - 1;
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

    pcout << getMultiBlockInfo(lattice) << std::endl;

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Initilization of rho and u." << endl;
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0 , u0);

    /*(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint nx, plint ny,
    Array<plint,3> position, plint radius, plint length, plint thickness)*/
    Array<plint,3> position(nx/2, ny/2, position_duct_z);
    plint thickness_duct = 2;
    build_duct(lattice, nx, ny, position, radius, length_duct, thickness_duct, omega);

    T rhoBar_target = 0;
    //const T mach_number = 0.2;
    //const T velocity_flow = mach_number*lattice_speed_sound;
    Array<T,3> j_target(0, 0, 0);
    T size_anechoic_buffer = 30;
    defineAnechoicMRTBoards(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;

    pcout << "Simulation begins" << endl;


    // --------------------------------------------------------
    // Setting probes
    // half size
    plint radius_probe = (radius - 1)/sqrt(2);
    // Two microphones 3r
    plint position_z_3r_p1 = position[2] + length_duct - 3*radius + 20;
    Box3D surface_probe_3r_p1(nx/2 - (radius_probe)/sqrt(2), 
            nx/2 + (radius_probe)/sqrt(2), 
            ny/2 - (radius_probe)/sqrt(2), 
            ny/2 + (radius_probe)/sqrt(2),
            position_z_3r_p1, position_z_3r_p1);

    plint position_z_3r_p2 = position_z_3r_p1 + 5;
    Box3D surface_probe_3r_p2(nx/2 - (radius_probe)/sqrt(2), 
            nx/2 + (radius_probe)/sqrt(2), 
            ny/2 - (radius_probe)/sqrt(2), 
            ny/2 + (radius_probe)/sqrt(2),
            position_z_3r_p2, position_z_3r_p2);
    
    // Two microphones 6r
    plint position_z_6r_p1 = position[2] + length_duct - 6*radius + 20;
    Box3D surface_probe_6r_p1(nx/2 - (radius_probe)/sqrt(2), 
            nx/2 + (radius_probe)/sqrt(2), 
            ny/2 - (radius_probe)/sqrt(2), 
            ny/2 + (radius_probe)/sqrt(2),
            position_z_6r_p1, position_z_6r_p1);

    plint position_z_6r_p2 = position_z_6r_p1 + 5;
    Box3D surface_probe_6r_p2(nx/2 - (radius_probe)/sqrt(2), 
            nx/2 + (radius_probe)/sqrt(2), 
            ny/2 - (radius_probe)/sqrt(2), 
            ny/2 + (radius_probe)/sqrt(2),
            position_z_6r_p2, position_z_6r_p2);

    // Boca of duct
    plint position_z_boca = position[2] + length_duct + 20;
    Box3D surface_probe_boca(nx/2 - (radius_probe)/sqrt(2), 
            nx/2 + (radius_probe)/sqrt(2), 
            ny/2 - (radius_probe)/sqrt(2), 
            ny/2 + (radius_probe)/sqrt(2),
            position_z_boca, position_z_boca);

    // Measurement in a boca
    std::string pressures_boca_string = fNameOut+"/history_pressures_boca.dat";
    char to_char_pressures_boca[1024];
    strcpy(to_char_pressures_boca, pressures_boca_string.c_str());
    plb_ofstream history_pressures_boca(to_char_pressures_boca);

    std::string velocities_boca_string = fNameOut+"/history_velocities_boca.dat";
    char to_char_velocities_boca[1024];
    strcpy(to_char_velocities_boca, velocities_boca_string.c_str());
    plb_ofstream history_velocities_boca(to_char_velocities_boca);

    // Measurement in 3r
    std::string pressures_3r_p1_string = fNameOut+"/history_pressures_3r_p1.dat";
    char to_char_pressures_3r_p1[1024];
    strcpy(to_char_pressures_3r_p1, pressures_3r_p1_string.c_str());
    plb_ofstream history_pressures_3r_p1(to_char_pressures_3r_p1);

    std::string pressures_3r_p2_string = fNameOut+"/history_pressures_3r_p2.dat";
    char to_char_pressures_3r_p2[1024];
    strcpy(to_char_pressures_3r_p2, pressures_3r_p2_string.c_str());
    plb_ofstream history_pressures_3r_p2(to_char_pressures_3r_p2);

    std::string velocities_3r_p1_string = fNameOut+"/history_velocities_3r_p1.dat";
    char to_char_velocities_3r_p1[1024];
    strcpy(to_char_velocities_3r_p1, velocities_3r_p1_string.c_str());
    plb_ofstream history_velocities_3r_p1(to_char_velocities_3r_p1);

    std::string velocities_3r_p2_string = fNameOut+"/history_velocities_3r_p2.dat";
    char to_char_velocities_3r_p2[1024];
    strcpy(to_char_velocities_3r_p2, velocities_3r_p2_string.c_str());
    plb_ofstream history_velocities_3r_p2(to_char_velocities_3r_p2);

    // Measurement in 6r
    std::string pressures_6r_p1_string = fNameOut+"/history_pressures_6r_p1.dat";
    char to_char_pressures_6r_p1[1024];
    strcpy(to_char_pressures_6r_p1, pressures_6r_p1_string.c_str());
    plb_ofstream history_pressures_6r_p1(to_char_pressures_6r_p1);

    std::string pressures_6r_p2_string = fNameOut+"/history_pressures_6r_p2.dat";
    char to_char_pressures_6r_p2[1024];
    strcpy(to_char_pressures_6r_p2, pressures_6r_p2_string.c_str());
    plb_ofstream history_pressures_6r_p2(to_char_pressures_6r_p2);

    std::string velocities_6r_p1_string = fNameOut+"/history_velocities_6r_p1.dat";
    char to_char_velocities_6r_p1[1024];
    strcpy(to_char_velocities_6r_p1, velocities_6r_p1_string.c_str());
    plb_ofstream history_velocities_6r_p1(to_char_velocities_6r_p1);

    std::string velocities_6r_p2_string = fNameOut+"/history_velocities_6r_p2.dat";
    char to_char_velocities_6r_p2[1024];
    strcpy(to_char_velocities_6r_p2, velocities_6r_p2_string.c_str());
    plb_ofstream history_velocities_6r_p2(to_char_velocities_6r_p2);

    // Informations about simulation    
    t = clock();
    std::string AllSimulationInfo_string = fNameOut + "/AllSimulationInfo.txt";
    char to_char_AllSimulationInfo[1024];
    strcpy(to_char_AllSimulationInfo, AllSimulationInfo_string.c_str());
    plb_ofstream AllSimulationInfo(to_char_AllSimulationInfo);

    std::string title = "\nSimulacao para testar a discretizacao. Sera que com mais elementos na na boca do duto dá certo?\n"; 
    AllSimulationInfo << endl
    << title << endl
    << "Dados da simulação" << endl
    << "Lattice:" << endl << endl 
    << "nx: " << nx << " ny: " << ny << " nz: " << nz << endl
    << " omega: " << omega << endl << endl
    << "Tempos: " << endl
    << "Total Time step: " << maxT << endl
    << "Discretizacao: " << radius/thickness_duct << endl
    << "Tamanho duto: " << length_duct << endl
    << "Posicao do duto: " << position[2] << endl;

    for (plint iT=0; iT<maxT; ++iT){
        if (iT <= maxT_final_source){
            
            plint total_signals = 20;
            T chirp_hand = get_linear_chirp_AZ(ka_max, total_signals, maxT_final_source, iT, drho, radius);
            
            Box3D place_source(position[0] - source_radius/sqrt(2), 
                position[0] + source_radius/sqrt(2), 
                position[1] - source_radius/sqrt(2), 
                position[1] + source_radius/sqrt(2), 
                position[2] + 2 + 19, position[2] + 2 + 20);
            //Box3D impulse(centerLB[0] + 10, centerLB[0] + 10, ny/2, ny/2, nz/2, nz/2);
            initializeAtEquilibrium(lattice, place_source, chirp_hand, u0);
        }else{
            Box3D place_source(position[0] - source_radius/sqrt(2), 
                position[0] + source_radius/sqrt(2), 
                position[1] - source_radius/sqrt(2), 
                position[1] + source_radius/sqrt(2), 
                position[2] + 2 + 19, position[2] + 2 + 20);
            initializeAtEquilibrium(lattice, place_source, rho0, u0);
        }


        if (iT % 10 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
        }

        if (iT % 1000 == 0) {
            //writeGifs(lattice,iT);
            writeVTK(lattice, iT);
        }

        // extract values of pressure and velocities
        // Extract from boca
        history_pressures_boca << setprecision(10) << (computeAverageDensity(lattice, surface_probe_boca) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_boca(plb::computeVelocityComponent(lattice, surface_probe_boca, 2));
        history_velocities_boca << setprecision(10) <<
        computeAverage(*velocity_boca, surface_probe_boca) << endl;
        // -----------------------------------------
        // Extract from 3r
        history_pressures_3r_p1 << setprecision(10) << (computeAverageDensity(lattice, surface_probe_3r_p1) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_3r_p1(plb::computeVelocityComponent(lattice, surface_probe_3r_p1, 2));
        history_velocities_3r_p1 << setprecision(10) <<
        computeAverage(*velocity_3r_p1, surface_probe_3r_p1) << endl;

        history_pressures_3r_p2 << setprecision(10) << (computeAverageDensity(lattice, surface_probe_3r_p2) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_3r_p2(plb::computeVelocityComponent(lattice, surface_probe_3r_p2, 2));
        history_velocities_3r_p2 << setprecision(10) <<
        computeAverage(*velocity_3r_p2, surface_probe_3r_p2) << endl;
        // -----------------------------------------
        // Extract from 6r
        history_pressures_6r_p1 << setprecision(10) << (computeAverageDensity(lattice, surface_probe_6r_p1) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_6r_p1(plb::computeVelocityComponent(lattice, surface_probe_6r_p1, 2));
        history_velocities_6r_p1 << setprecision(10) <<
        computeAverage(*velocity_6r_p1, surface_probe_6r_p1) << endl;

        history_pressures_6r_p2 << setprecision(10) << (computeAverageDensity(lattice, surface_probe_6r_p2) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_6r_p2(plb::computeVelocityComponent(lattice, surface_probe_6r_p2, 2));
        history_velocities_6r_p2 << setprecision(10) <<
        computeAverage(*velocity_6r_p2, surface_probe_6r_p2) << endl;
        // -----------------------------------------

        lattice.collideAndStream();
    }

    t = (clock() - t)/CLOCKS_PER_SEC;
    AllSimulationInfo << endl << "Execution time: " << t << " segundos" << endl;

    pcout << "End of simulation at iteration " << endl;
}