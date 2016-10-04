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

void build_duct(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint nx, plint ny,
    Array<plint,3> position, plint radius, plint length, plint thickness, T omega){
    length += 4;
    // Duct is constructed along the Z direction
    //plint size_square = 50;
    plint size_square = 2*radius;
    plint radius_intern = radius - thickness;
    for (plint x = position[0] - radius; x < nx/2 + size_square/2; ++x){
        for (plint y = position[1] - radius; y < ny/2 + size_square/2; ++y){
            for (plint z = position[2]; z < length + position[2]; ++z){

                if (radius*radius > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2)){
                    //pcout << "passou" << endl;
                    DotList3D points_to_aplly_dynamics;
                    points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                    defineDynamics(lattice, points_to_aplly_dynamics, new BounceBack<T,DESCRIPTOR>(0));
                }
                // extrude
                if (radius_intern*radius_intern > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2) && z > position[2] + 2){
                    //pcout << "passou" << endl;
                    DotList3D points_to_aplly_dynamics;
                    points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                    defineDynamics(lattice, points_to_aplly_dynamics, new BackgroundDynamics(omega));
                }

            }
        }
    }

}

T computeMeanVelocityComponent(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Box3D surface_probe, plint component){
    Array<T,6> positions = surface_probe.to_plbArray();
    T mean_velocity = 0;
    

    for (plint x = positions[0]; x <= positions[1]; ++x){
        for (plint y = positions[2]; y <= positions[3]; ++y){
            for (plint z = positions[4]; z <= positions[5]; ++z){
                Array<T,3> velocities;
                lattice.get(x, y, z).computeVelocity(velocities);
                mean_velocity += velocities[component];
                //pcout << "entrei" <<  endl;
                //pcout << mean_velocity <<  endl;
            }
        }        
    }
    mean_velocity = mean_velocity/surface_probe.nCells();    
    
    return mean_velocity;
}

int main(int argc, char **argv){
    plbInit(&argc, &argv);

    //const plint length_domain = 420;
    const plint radius = 20;
    const plint diameter = 2*radius;
    //const plint length_domain = 150;
    const plint nx = 6*diameter;
    const plint ny = 6*diameter;
    const plint position_duct_z = 30;
    const plint nz = 8*diameter;
    const T lattice_speed_sound = 1/sqrt(3);
    const T omega = 1.985;
    const plint maxT = 20000;

    //const plint maxT = 2*120/lattice_speed_sound;
    const plint maxT_final_source = maxT - nz*sqrt(3);
    const T ka_max = 2.5;
    const T ka_min = 0;
    const T cs2 = lattice_speed_sound*lattice_speed_sound;
    clock_t t;

    Array<T,3> u0(0, 0, 0);

    std::string fNameOut = "tmp";
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

    /*plint size_square = 10;
    Box3D square(
    35, 35 + size_square,
    ny/2 - size_square/2, ny/2 + size_square/2, 
    nz/2 - size_square/2, nz/2 + size_square/2);
    defineDynamics(lattice, square, new BounceBack<T,DESCRIPTOR>((T)0));*/

    /*(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint nx, plint ny,
    Array<plint,3> position, plint radius, plint length, plint thickness)*/
    Array<plint,3> position(nx/2, ny/2, position_duct_z);
    plint length_duct = 120;
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

    // Setting probes
    plint radius_probe = (radius - 1)/sqrt(2);
    // half size
    plint position_z_3r = position[2] + length_duct - 3*radius;
    Box3D surface_probe_3r(nx/2 - radius_probe/sqrt(2), 
            nx/2 + radius_probe/sqrt(2), 
            ny/2 - radius_probe/sqrt(2), 
            ny/2 + radius_probe/sqrt(2),
            position_z_3r, position_z_3r);

    //plint position_z_4r = position[2] + length_duct - 4*radius;
    plint position_z_4r = position[2] + length_duct;
    Box3D surface_probe_4r(nx/2 - radius_probe/sqrt(2), 
            nx/2 + radius_probe/sqrt(2), 
            ny/2 - radius_probe/sqrt(2), 
            ny/2 + radius_probe/sqrt(2),
            position_z_4r, position_z_4r);

    plint position_z_6r = position[2] + length_duct - 6*radius;
    Box3D surface_probe_6r(nx/2 - radius_probe/sqrt(2), 
            nx/2 + radius_probe/sqrt(2), 
            ny/2 - radius_probe/sqrt(2), 
            ny/2 + radius_probe/sqrt(2),
            position_z_6r, position_z_6r);
    plb_ofstream history_pressures_3r("tmp/history_pressures_3r.dat");
    plb_ofstream history_pressures_4r("tmp/history_pressures_boca.dat");
    plb_ofstream history_pressures_6r("tmp/history_pressures_6r.dat");
    plb_ofstream history_velocities_3r("tmp/history_velocities_3r.dat");
    plb_ofstream history_velocities_4r("tmp/history_velocities_boca.dat");
    plb_ofstream history_velocities_6r("tmp/history_velocities_6r.dat");

    plb_ofstream history_pressures_3r_point("tmp/history_pressures_3r_point.dat");
    plb_ofstream history_pressures_4r_point("tmp/history_pressures_boca_point.dat");
    plb_ofstream history_pressures_6r_point("tmp/history_pressures_6r_point.dat");
    plb_ofstream history_velocities_3r_point("tmp/history_velocities_3r_point.dat");
    plb_ofstream history_velocities_4r_point("tmp/history_velocities_boca_point.dat");
    plb_ofstream history_velocities_6r_point("tmp/history_velocities_6r_point.dat");

    t = clock();
    plb_ofstream AllSimulationInfo("./tmp/AllSimulationInfo.txt");

    AllSimulationInfo << endl 
    << "Dados da simulação" << endl
    << "Lattice:" << endl << endl 
    << "nx: " << nx << " ny: " << ny << " nz: " << nz << endl
    << " omega: " << omega << endl << endl
    << "Tempos: " << endl
    << "Total Time step: " << maxT << endl
    << "Discretizacao: " << radius/thickness_duct << endl
    << "Tamanho duto: " << length_duct << endl
    << "Posicao do duto: " << position[2] << endl
    ;
    for (plint iT=0; iT<maxT; ++iT){
        if (iT <= maxT_final_source){
            //drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            //drho*cos((lattice_speed_sound/radius)*(ka_max*((maxT-iT)/maxT)));
            
            T initial_frequency = ka_min*lattice_speed_sound/(2*M_PI*radius);
            T frequency_max_lattice = ka_max*lattice_speed_sound/(2*M_PI*radius);
            T variation_frequency = (frequency_max_lattice - initial_frequency)/maxT_final_source;
            T frequency_function = initial_frequency*iT + (variation_frequency*iT*iT)/2;
            T phase = 2*M_PI*frequency_function;
            T chirp_hand = 1. + drho*sin(phase);

            //T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            //Box3D impulse(nx/2, nx/2, ny/2, ny/2, position_z_3r, position_z_3r);
            plint source_radius = radius - 1;
            Box3D place_source(position[0] - source_radius/sqrt(2), 
                position[0] + source_radius/sqrt(2), 
                position[1] - source_radius/sqrt(2), 
                position[1] + source_radius/sqrt(2), 
                position[2] + 3, position[2] + 8);
            //Box3D impulse(centerLB[0] + 10, centerLB[0] + 10, ny/2, ny/2, nz/2, nz/2);
            initializeAtEquilibrium(lattice, place_source, chirp_hand, u0);
        }else{
            plint source_radius = radius - 1;
            Box3D place_source(position[0] - source_radius/sqrt(2), 
                position[0] + source_radius/sqrt(2), 
                position[1] - source_radius/sqrt(2), 
                position[1] + source_radius/sqrt(2), 
                position[2] + 3, position[2] + 8);
            initializeAtEquilibrium(lattice, place_source, rho0, u0);
        }

        if (iT % 100 == 0 && iT>0) {
            pcout << "Iteration " << iT << endl;
        }

        if (iT % 100 == 0) {
            //pcout << "Iteration " << iT << endl;
            //writeGifs(lattice,iT);
            writeVTK(lattice, iT);
        }

        // extract values of pressure and velocities
        history_pressures_3r << setprecision(10) << (computeAverageDensity(lattice, surface_probe_3r) - rho0)*cs2 << endl;
        history_pressures_4r << setprecision(10) << (computeAverageDensity(lattice, surface_probe_4r) - rho0)*cs2 << endl;
        history_pressures_6r << setprecision(10) << (computeAverageDensity(lattice, surface_probe_6r) - rho0)*cs2 << endl;


        history_velocities_3r << setprecision(10) << 
        computeMeanVelocityComponent(lattice, surface_probe_3r, 2)/lattice_speed_sound << endl;

        history_velocities_4r << setprecision(10) << 
        computeMeanVelocityComponent(lattice, surface_probe_4r, 2)/lattice_speed_sound << endl;

        history_velocities_6r << setprecision(10) << 
        computeMeanVelocityComponent(lattice, surface_probe_6r, 2)/lattice_speed_sound << endl;


        // extract values of pressure and velocities point
        history_pressures_3r_point << setprecision(10) << (lattice.get(nx/2, ny/2, position_z_3r).computeDensity() - rho0)*cs2 << endl;
        history_pressures_4r_point << setprecision(10) << (lattice.get(nx/2, ny/2, position_z_4r).computeDensity() - rho0)*cs2 << endl;
        history_pressures_6r_point << setprecision(10) << (lattice.get(nx/2, ny/2, position_z_6r).computeDensity() - rho0)*cs2 << endl;

        Array<T,3> velocities;
        lattice.get(nx/2, ny/2, position_z_3r).computeVelocity(velocities);
        history_velocities_3r_point << setprecision(10) << velocities[2]/lattice_speed_sound << endl;
        lattice.get(nx/2, ny/2, position_z_4r).computeVelocity(velocities);
        history_velocities_4r_point << setprecision(10) << velocities[2]/lattice_speed_sound << endl;
        lattice.get(nx/2, ny/2, position_z_6r).computeVelocity(velocities);
        history_velocities_6r_point << setprecision(10) << velocities[2]/lattice_speed_sound << endl;


        lattice.collideAndStream();

    }

     t = (clock() - t)/CLOCKS_PER_SEC;
    AllSimulationInfo << endl << "Execution time: " << t/60 << " minutos" << endl << endl;

    pcout << "End of simulation at iteration " << endl;
}
