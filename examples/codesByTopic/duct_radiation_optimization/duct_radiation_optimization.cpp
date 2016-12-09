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
const T drho = rho0/10;

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
    //T mach_number = 0.18;
    T mach_number = 0;
    T lattice_speed_sound = 1/sqrt(3);
    T velocity_flow = -mach_number*lattice_speed_sound;
    plint anechoic_size = 29;
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
                        Array<T,3> u0(0, 0, velocity_flow);
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

T get_linear_chirp(T ka_min, T ka_max, plint maxT_final_source, plint iT, T drho, T radius){
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
    total_signals = 2*total_signals;
    for (plint n_signal = 1; n_signal <= total_signals; n_signal++){
        T interval = ka_max/total_signals;
        T phase = (n_signal*interval*cs*iT)/(radius);
        chirp_value += drho*sin(phase);
    }
    return chirp_value;
}


void set_source(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Array<plint,3> position, 
    T chirp_hand, Array<T,3> u0, plint radius, plint radius_intern, plint nx, plint ny){

    plint size_square = 2*radius;
    Box3D impulse_local;
    for (plint x = position[0] - radius; x < nx/2 + size_square/2; ++x){
        for (plint y = position[1] - radius; y < ny/2 + size_square/2; ++y){
            // extrude
            if (radius_intern*radius_intern > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2)){
                Array<plint, 6> local_source_2(x, x, y, y, position[2] + 29, position[2] + 30);
                impulse_local.from_plbArray(local_source_2);
                Array<T,3> u_chirp_hand(0, 0, (chirp_hand-1)/(1/sqrt(3)));
                initializeAtEquilibrium(lattice, impulse_local, chirp_hand, u_chirp_hand);
            }
        }
    }
}

T compute_avarage_density(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Array<plint,3> position, 
    plint radius, plint radius_intern, plint nx, plint ny, plint position_z){

    plint number_cells = 0;
    T average_density = 0;
    plint size_square = 2*radius;
    Box3D local_box;
    for (plint x = position[0] - radius; x < nx/2 + size_square/2; ++x){
        for (plint y = position[1] - radius; y < ny/2 + size_square/2; ++y){
            // extrude
            if (radius_intern*radius_intern > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2)){
                Array<plint, 6> local(x, x, y, y, position_z, position_z);
                local_box.from_plbArray(local);
                average_density += computeAverageDensity(lattice, local_box);
                number_cells += 1;
            }
        }
    }

    average_density = average_density/number_cells;

    return average_density;
}

T compute_avarage_velocity(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Array<plint,3> position, 
    plint radius, plint radius_intern, plint nx, plint ny, plint position_z){

    plint number_cells = 0;
    T average_velocity = 0;
    plint size_square = 2*radius;
    Box3D local_box;
    for (plint x = position[0] - radius; x < nx/2 + size_square/2; ++x){
        for (plint y = position[1] - radius; y < ny/2 + size_square/2; ++y){
            // extrude
            if (radius_intern*radius_intern > (x-nx/2)*(x-nx/2) + (y-ny/2)*(y-ny/2)){
                Array<plint, 6> local(x, x, y, y, position_z, position_z);
                local_box.from_plbArray(local);
                std::auto_ptr<MultiScalarField3D<T> > velocity(plb::computeVelocityComponent(lattice, local_box, 2));
                average_velocity += computeAverage(*velocity, local_box);
                number_cells += 1;
            }
        }
    }

    average_velocity = average_velocity/number_cells;

    return average_velocity;
}


void set_nodynamics(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint nx, plint ny, plint off_set_z){
    Box3D place_nodynamics(0, nx - 1, 0, ny - 1, 0, off_set_z);
    defineDynamics(lattice, place_nodynamics, new NoDynamics<T,DESCRIPTOR>(0));
}

class Probe{
    private:
        Box3D location;
        string name_probe;
        plb_ofstream file_pressures;
        plb_ofstream file_velocities_x;
        plb_ofstream file_velocities_y;
        plb_ofstream file_velocities_z;
    public:
        Probe() : location(), name_probe(), file_pressures(), file_velocities_x(),
        file_velocities_y(), file_velocities_z(){ }
        Probe(Box3D location, string directory, string name_probe){

            directory = directory + "/" + name_probe;
            std::string command = "mkdir -p " + directory;
            char to_char_command[1024];
            strcpy(to_char_command, command.c_str());
            system(to_char_command);

            string pressures_string = directory + "/history_pressures_" + name_probe + ".dat";
            string velocities_x_string = directory + "/history_velocities_x_" + name_probe + ".dat";
            string velocities_y_string = directory + "/history_velocities_y_" + name_probe + ".dat";
            string velocities_z_string = directory + "/history_velocities_z_" + name_probe + ".dat";
            char to_char_pressures[1024];
            char to_char_velocities_x[1024];
            char to_char_velocities_y[1024];
            char to_char_velocities_z[1024];
            strcpy(to_char_pressures, pressures_string.c_str());
            strcpy(to_char_velocities_x, velocities_x_string.c_str());
            strcpy(to_char_velocities_y, velocities_y_string.c_str());
            strcpy(to_char_velocities_z, velocities_z_string.c_str());

            this->file_pressures.open(to_char_pressures);
            this->file_velocities_x.open(to_char_velocities_x);
            this->file_velocities_y.open(to_char_velocities_y);
            this->file_velocities_z.open(to_char_velocities_z);
            this->location = location;
            this->name_probe = name_probe;
        }

        string get_name_probe(){
            return this->name_probe;
        }

        void set_properties(Box3D location, string directory, string name_probe){
            directory = directory + "/" + name_probe;
            std::string command = "mkdir -p " + directory;
            char to_char_command[1024];
            strcpy(to_char_command, command.c_str());
            system(to_char_command);

            string pressures_string = directory + "/history_pressures_" + name_probe + ".dat";
            string velocities_x_string = directory + "/history_velocities_x_" + name_probe + ".dat";
            string velocities_y_string = directory + "/history_velocities_y_" + name_probe + ".dat";
            string velocities_z_string = directory + "/history_velocities_z_" + name_probe + ".dat";
            char to_char_pressures[1024];
            char to_char_velocities_x[1024];
            char to_char_velocities_y[1024];
            char to_char_velocities_z[1024];
            strcpy(to_char_pressures, pressures_string.c_str());
            strcpy(to_char_velocities_x, velocities_x_string.c_str());
            strcpy(to_char_velocities_y, velocities_y_string.c_str());
            strcpy(to_char_velocities_z, velocities_z_string.c_str());

            this->file_pressures.open(to_char_pressures);
            this->file_velocities_x.open(to_char_velocities_x);
            this->file_velocities_y.open(to_char_velocities_y);
            this->file_velocities_z.open(to_char_velocities_z);
            this->location = location;
            this->name_probe = name_probe;   
        }

        void save_point(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T rho0, T cs2){
            file_pressures << setprecision(10) << (computeAverageDensity(lattice, this->location) - rho0)*cs2 << endl;
            std::auto_ptr<MultiScalarField3D<T> > velocity_x(plb::computeVelocityComponent(lattice, this->location, 0));
            file_velocities_x << setprecision(10) << computeAverage(*velocity_x, this->location) << endl;
            std::auto_ptr<MultiScalarField3D<T> > velocity_y(plb::computeVelocityComponent(lattice, this->location, 1));
            file_velocities_y << setprecision(10) << computeAverage(*velocity_y, this->location) << endl;
            std::auto_ptr<MultiScalarField3D<T> > velocity_z(plb::computeVelocityComponent(lattice, this->location, 2));
            file_velocities_z << setprecision(10) << computeAverage(*velocity_z, this->location) << endl;
        }
};

class Two_Microphones{
    private:
        Probe microphone_1;
        Probe microphone_2;
    public:
        Two_Microphones(plint radius, plint microphone_distance, 
            plint length_duct, Array<plint,3> position_duct, 
            std::string fNameOut, std::string name, plint distance, 
            plint nx, plint ny, plint nz){

            plint radius_probe = (radius - 1)/sqrt(2);
            plint position_z = position_duct[2] + length_duct - distance;
            
            Box3D surface_probe_p1(nx/2 - (radius_probe)/sqrt(2), 
                    nx/2 + (radius_probe)/sqrt(2), 
                    ny/2 - (radius_probe)/sqrt(2), 
                    ny/2 + (radius_probe)/sqrt(2),
                    position_z, position_z);
            std::string name_p1 = name + "_p1";
            Probe probe_p1(surface_probe_p1, fNameOut, name_p1);
            this->microphone_1.set_properties(surface_probe_p1, fNameOut, name_p1);

            Box3D surface_probe_p2(nx/2 - (radius_probe)/sqrt(2), 
                    nx/2 + (radius_probe)/sqrt(2), 
                    ny/2 - (radius_probe)/sqrt(2), 
                    ny/2 + (radius_probe)/sqrt(2),
                    position_z + microphone_distance,
                    position_z + microphone_distance);
            std::string name_p2 = name + "_p2";
            this->microphone_2.set_properties(surface_probe_p2, fNameOut, name_p2);
        }

        void save_point(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T rho0, T cs2){
            this->microphone_1.save_point(lattice, rho0, cs2);
            this->microphone_2.save_point(lattice, rho0, cs2);   
        }
};

class System_Abom_Measurement{
    private:
        std::vector<Box3D> microphones_positions; 
        plb_ofstream file_pressures;
        plb_ofstream file_velocities_x;
        plb_ofstream file_velocities_y;
        plb_ofstream file_velocities_z;
    public:
        System_Abom_Measurement(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Array<plint,3> position_duct, 
            plint radius, plint distance_group_A, plint distance_group_B, string directory){

            string name_probe = "system_abom_measurement_data";
            directory = directory + "/" + name_probe;
            std::string command = "mkdir -p " + directory;
            char to_char_command[1024];
            strcpy(to_char_command, command.c_str());
            system(to_char_command);

            // normaly distance_group_A = 28
            // normaly distance_group_B = 135

            const plint nx = lattice.getNx();
            const plint ny = lattice.getNy();
            const plint nz = lattice.getNz();
            const plint diameter = 2*radius;

            string pressures_string = directory + "/history_pressures_" + name_probe + ".dat";
            string velocities_x_string = directory + "/history_velocities_x_" + name_probe + ".dat";
            string velocities_y_string = directory + "/history_velocities_y_" + name_probe + ".dat";
            string velocities_z_string = directory + "/history_velocities_z_" + name_probe + ".dat";
            char to_char_pressures[1024];
            char to_char_velocities_x[1024];
            char to_char_velocities_y[1024];
            char to_char_velocities_z[1024];
            strcpy(to_char_pressures, pressures_string.c_str());
            strcpy(to_char_velocities_x, velocities_x_string.c_str());
            strcpy(to_char_velocities_y, velocities_y_string.c_str());
            strcpy(to_char_velocities_z, velocities_z_string.c_str());

            this->file_pressures.open(to_char_pressures);
            this->file_velocities_x.open(to_char_velocities_x);
            this->file_velocities_y.open(to_char_velocities_y);
            this->file_velocities_z.open(to_char_velocities_z);

            // Group A
            plint position_microphone_1 = position_duct[2] + 30 + 3*diameter;
            Box3D microphone_1_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_1, position_microphone_1);
            microphones_positions.push_back(microphone_1_position);

            plint position_microphone_2 = position_microphone_1 + distance_group_A;
            Box3D microphone_2_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_2, position_microphone_2);
            microphones_positions.push_back(microphone_2_position);

            plint position_microphone_3 = position_microphone_2 + distance_group_A;
            Box3D microphone_3_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_3, position_microphone_3);
            microphones_positions.push_back(microphone_3_position);

            plint position_microphone_4 = position_microphone_3 + distance_group_A;
            Box3D microphone_4_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_4, position_microphone_4);
            microphones_positions.push_back(microphone_4_position);

            plint position_microphone_5 = position_microphone_4 + distance_group_A;
            Box3D microphone_5_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_5, position_microphone_5);
            microphones_positions.push_back(microphone_5_position);

            plint position_microphone_6 = position_microphone_5 + distance_group_A;
            Box3D microphone_6_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_6, position_microphone_6);
            microphones_positions.push_back(microphone_6_position);

            // Group B
            plint position_microphone_7 = position_duct[2] + 30 + 3*diameter;
            Box3D microphone_7_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_7, position_microphone_7);
            microphones_positions.push_back(microphone_7_position);

            plint position_microphone_8 = position_microphone_7 + distance_group_B;
            Box3D microphone_8_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_8, position_microphone_8);
            microphones_positions.push_back(microphone_8_position);

            plint position_microphone_9 = position_microphone_8 + distance_group_B;
            Box3D microphone_9_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_9, position_microphone_9);
            microphones_positions.push_back(microphone_9_position);

            plint position_microphone_10 = position_microphone_9 + distance_group_B;
            Box3D microphone_10_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_10, position_microphone_10);
            microphones_positions.push_back(microphone_10_position);

            plint position_microphone_11 = position_microphone_10 + distance_group_B;
            Box3D microphone_11_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_11, position_microphone_11);
            microphones_positions.push_back(microphone_11_position);

            plint position_microphone_12 = position_microphone_11 + distance_group_B;
            Box3D microphone_12_position(position_duct[0], position_duct[0],
                position_duct[1], position_duct[1],
                position_microphone_12, position_microphone_12);
            microphones_positions.push_back(microphone_12_position);

            pcout << "MICROPHONES" << endl;
            for (int i = 0; i < 12; ++i)
            {
                Box3D to_see = microphones_positions[i];
                Array<plint, 6> test = to_see.to_plbArray();
                pcout << test[5] << endl;
            }
        }

    void save_point(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T rho0, T cs2){

        for (int mic = 0; mic < 11; mic++){
            file_pressures << setprecision(10) << (computeAverageDensity(lattice, this->microphones_positions[mic]) - rho0)*cs2 << " ";
            std::auto_ptr<MultiScalarField3D<T> > velocity_x(plb::computeVelocityComponent(lattice, this->microphones_positions[mic], 0));
            file_velocities_x << setprecision(10) << computeAverage(*velocity_x, this->microphones_positions[mic]) << " ";
            std::auto_ptr<MultiScalarField3D<T> > velocity_y(plb::computeVelocityComponent(lattice, this->microphones_positions[mic], 1));
            file_velocities_y << setprecision(10) << computeAverage(*velocity_y, this->microphones_positions[mic]) << " ";
            std::auto_ptr<MultiScalarField3D<T> > velocity_z(plb::computeVelocityComponent(lattice, this->microphones_positions[mic], 2));
            file_velocities_z << setprecision(10) << computeAverage(*velocity_z, this->microphones_positions[mic]) << " ";
        }

        file_pressures << setprecision(10) << (computeAverageDensity(lattice, this->microphones_positions[11]) - rho0)*cs2 << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_x(plb::computeVelocityComponent(lattice, this->microphones_positions[11], 0));
        file_velocities_x << setprecision(10) << computeAverage(*velocity_x, this->microphones_positions[11]) << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_y(plb::computeVelocityComponent(lattice, this->microphones_positions[11], 1));
        file_velocities_y << setprecision(10) << computeAverage(*velocity_y, this->microphones_positions[11]) << endl;
        std::auto_ptr<MultiScalarField3D<T> > velocity_z(plb::computeVelocityComponent(lattice, this->microphones_positions[11], 2));
        file_velocities_z << setprecision(10) << computeAverage(*velocity_z, this->microphones_positions[11]) << endl;
    }
};

int main(int argc, char **argv){
    plbInit(&argc, &argv);

    const plint radius = 20;
    const plint diameter = 2*radius;
    const plint nx = 6*diameter + 60;
    const plint ny = 6*diameter + 60;
    const plint position_duct_z = 0;
    //const plint length_duct = 10*diameter + 30;
    const plint length_duct = 30 + 120 + 5*113 + 3*diameter;
    const plint nz = length_duct + 3*diameter + 30;
    const T lattice_speed_sound = 1/sqrt(3);
    const T omega = 1.985;
    const plint maxT = pow(2,13) + nz*sqrt(3);
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
    //const T mach_number = 0.2;
    //const T velocity_flow = mach_number*lattice_speed_sound;
    Array<T,3> j_target(0, 0, 0);
    T size_anechoic_buffer = 30;
    defineAnechoicMRTBoards_limited(nx, ny, nz, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target, j_target, j_target,
      rhoBar_target, off_set_z);

    build_duct(lattice, nx, ny, position, radius, length_duct, thickness_duct, omega);

    lattice.initialize();

    pcout << std::endl << "Voxelizing the domain." << std::endl;

    pcout << "Simulation begins" << endl;

    // Setting probes ------------------------------------------
    plint distance_group_A = 28;
    plint distance_group_B = 113;
    System_Abom_Measurement system_abom_measurement(lattice, position, radius, 
        distance_group_A, distance_group_B, fNameOut);
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
    
    std::string title = "\nTENTANDO DE NOVO A TECNICA DO ABOM.\n"; 
    
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
        if (iT <= maxT_final_source){
            plint total_signals = 20;
            T chirp_hand = get_linear_chirp(ka_min, ka_max, maxT_final_source, iT, drho, radius);
            //T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            history_signal_in << setprecision(10) << chirp_hand << endl;
            set_source(lattice, position, chirp_hand, u0, radius, radius_intern, nx, ny);
        }else{
            set_source(lattice, position, rho0, u0, radius, radius_intern, nx, ny);
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

        if (iT % 50 == 0) {
            //writeGifs(lattice,iT);
            //writeVTK(lattice, iT);
        }

        // extract values of pressure and velocities
        system_abom_measurement.save_point(lattice, rho0, cs2);

        lattice.collideAndStream();
    }

    t = (clock() - t)/CLOCKS_PER_SEC;
    AllSimulationInfo << endl << "Execution time: " << t << " segundos" << endl;

    pcout << "End of simulation at iteration " << endl;
}