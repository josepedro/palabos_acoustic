/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Code 1.1 in the Palabos tutorial
 */

#include "palabos2D.h"
#ifndef PLB_PRECOMPILED // Unless precompiled version is used,
#include "palabos2D.hh"   // include full template code
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <string> 

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
//#define DESCRIPTOR plb::descriptors::D2Q9Descriptor
#define DESCRIPTOR MRTD2Q9Descriptor
typedef MRTdynamics<T,DESCRIPTOR> BackgroundDynamics;
typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;

// ---------------------------------------------
// Includes of acoustics resources
#include "acoustics/acoustics2D.h"
using namespace plb_acoustics_2D;
// ---------------------------------------------

const T c_o = 343;
const T density_phys = 1.22;
const T freq_phys = 500; // Hz
//const T delta_x = 2*10e-5;
const T delta_x = 2*10e-4;
const T NPS = 100;
const T rho0 = 1.;
const Array<T,2> u0((T) 0,(T)0);
const T delta_t = 5.66*10e-9;
const T As_lattice = (2*10e-3)/delta_x;
const T d_thickness_phys = 10e-3;
const plint d_thickness_lattice = d_thickness_phys/delta_x;
const T size_anechoic_buffer = 50;
const plint HH = (120*10e-3)/delta_x;
const plint H = (150*10e-3)/delta_x;
const plint ny = 2*size_anechoic_buffer + HH + d_thickness_lattice + H;
const plint At = (35*10e-3)/delta_x;
const plint Af = ((120*10e-3)/delta_x - At)/2;
const plint nx = 2*size_anechoic_buffer + 2*Af + At;       // Choice of lattice dimensions.
const T cs = 1/sqrt(3);
const T cs2 = cs*cs;
const T zeta = c_o/cs;
const plint period_cicle_duration = HH/cs;
const plint transient_time = 20*period_cicle_duration;
const plint record_time_points = 5*period_cicle_duration;
const plint maxIter = 2*(transient_time + record_time_points + 10000);
const T omega = 1.97;
const T tau = 1/omega;
const T kinematic_viscosity_lattice = cs2*(tau - 0.5);
const T kinematic_viscosity_phys = zeta*delta_x*kinematic_viscosity_lattice;
const T dynamic_viscosity_phys = kinematic_viscosity_phys/density_phys;
const T reynolds_number = (density_phys*c_o*d_thickness_phys)/dynamic_viscosity_phys;
const T mach_number = 0;
const T velocity_flow = mach_number*cs;

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    std::string fNameOut = currentDateTime() + "+tmp";
    std::string command = "mkdir -p " + fNameOut;
    char to_char_command[1024];
    strcpy(to_char_command, command.c_str());
    system(to_char_command);

    global::directories().setOutputDir(fNameOut+"/");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new AnechoicBackgroundDynamics(omega));
    defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;
    pcout << "omega: " << omega << std::endl;
    pcout << "velocity_flow: " << velocity_flow << std::endl;
    pcout << "Total iteration: " << maxIter << std::endl;
    pcout << getMultiBlockInfo(lattice) << std::endl;
    pcout << "Total memory RAM required (Gb): " << T (((nx*ny*19*8)/1000)/1000)/1000 << std::endl;

    std::string file_string_info = fNameOut + "/reynolds_dynamic_viscosity_phys.dat";
    char to_char_file_string_info[1024];
    strcpy(to_char_file_string_info, file_string_info.c_str());
    plb_ofstream info_reynols_dynamic_viscosity_phys(to_char_file_string_info);
    info_reynols_dynamic_viscosity_phys << "Reynolds Number: "  <<  setprecision(10) << reynolds_number  
    << endl << "Dynamic Viscosity Physical: " << setprecision(10) << dynamic_viscosity_phys << endl;

    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    
    lattice.initialize();

    // Anechoic Condition
    T rhoBar_target = 0;
    Array<T,2> j_target(0./std::sqrt(3), 0.0/std::sqrt(3));
    // Define Anechoic Boards
    defineAnechoicMRTBoards(nx, ny, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target,
      rhoBar_target, rhoBar_target, rhoBar_target, rhoBar_target);

    Box2D square_1(0, Af + size_anechoic_buffer, 0, H + size_anechoic_buffer + d_thickness_lattice);
    defineDynamics(lattice, square_1, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_2(0, nx - 1, 0, size_anechoic_buffer - 1);
    defineDynamics(lattice, square_2, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_3(nx - Af - size_anechoic_buffer, nx - 1, 0,  H + size_anechoic_buffer + d_thickness_lattice);
    defineDynamics(lattice, square_3, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_4(0, size_anechoic_buffer + Af + (At - As_lattice)/2,
    size_anechoic_buffer + H, size_anechoic_buffer + H  + d_thickness_lattice);
    defineDynamics(lattice, square_4, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_5(size_anechoic_buffer + Af + (At - As_lattice)/2 + As_lattice, nx - 1, 
    size_anechoic_buffer + H, size_anechoic_buffer + H  + d_thickness_lattice);
    defineDynamics(lattice, square_5, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D above(0, nx-1, ny - size_anechoic_buffer - 1, ny - 1);
    defineDynamics(lattice, above, new BounceBack<T,DESCRIPTOR>((T)0));

    // Main loop over time iterations.
    /*Box2D signal_input_place(size_anechoic_buffer + Af, nx - Af - size_anechoic_buffer,
    ny - size_anechoic_buffer - 2, ny  - size_anechoic_buffer - 2);*/
    Box2D signal_input_place(nx - size_anechoic_buffer - 1, nx - size_anechoic_buffer - 1,
    size_anechoic_buffer + H  + d_thickness_lattice, ny  - size_anechoic_buffer - 2);
    //Box2D local_to_extract(Af + size_anechoic_buffer, Af + size_anechoic_buffer + At,
    //size_anechoic_buffer, 2*H+d_thickness_lattice + size_anechoic_buffer);
    Box2D local_to_extract(Af + size_anechoic_buffer, Af + size_anechoic_buffer + At,
    size_anechoic_buffer + H - 5*d_thickness_lattice, size_anechoic_buffer + H  + 6*d_thickness_lattice);

    info_reynols_dynamic_viscosity_phys << "Nx local: "  <<  setprecision(10) << local_to_extract.getNx()  
    << endl << "Ny local: " << setprecision(10) << local_to_extract.getNy()
    << endl << "Total number iterations: " << maxIter + 1 << endl;

    Box2D pressures_1_extract_local(size_anechoic_buffer + Af - (15*10e-3)/delta_x,size_anechoic_buffer + Af - (15*10e-3)/delta_x,
    size_anechoic_buffer + H  + d_thickness_lattice + 2,size_anechoic_buffer + H  + d_thickness_lattice + 2);
    std::string file_string_pressures_1 = fNameOut + "/pressures_1.dat";
    char to_char_file_string_pressures_1[1024];
    strcpy(to_char_file_string_pressures_1, file_string_pressures_1.c_str());
    plb_ofstream pressures_1(to_char_file_string_pressures_1);

    Box2D pressures_2_extract_local(size_anechoic_buffer + Af + At/2,size_anechoic_buffer + Af + At/2, 
        size_anechoic_buffer + 2, size_anechoic_buffer + 2);
    std::string file_string_pressures_2 = fNameOut + "/pressures_2.dat";
    char to_char_file_string_pressures_2[1024];
    strcpy(to_char_file_string_pressures_2, file_string_pressures_2.c_str());
    plb_ofstream pressures_2(to_char_file_string_pressures_2);
    for (plint iT = 0; iT <= maxIter; iT++){

        // excitation signal
        if (iT >= 0){
            T rho_changing = get_tone_pure_omega(freq_phys, iT, NPS, delta_x);
            initializeAtEquilibrium(lattice, signal_input_place, rho_changing, u0);
        }

        if (iT%50==0) { 
            pcout << "iT= " << iT << endl;

            

            pcout << " energy ="
            << setprecision(10) << getStoredAverageEnergy<T>(lattice)
            << " rho ="
            << getStoredAverageDensity<T>(lattice)
            << " max_velocity ="
            << setprecision(10) << (getStoredMaxVelocity<T>(lattice))///cs
            << endl;
        }

	if(iT%100==0){
		//writeVTK(lattice, iT, local_to_extract);
        writeVTK(lattice, iT, lattice.getBoundingBox());
	}
       
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
        
        pressures_1 << setprecision(10) << (computeAverageDensity(lattice, pressures_1_extract_local) - rho0)*cs2 << endl;
        pressures_2 << setprecision(10) << (computeAverageDensity(lattice, pressures_2_extract_local) - rho0)*cs2 << endl;
    }

}
