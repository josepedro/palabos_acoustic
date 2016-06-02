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

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

// ---------------------------------------------
// Includes of acoustics resources
#include "acoustics/acoustics2D.h"
using namespace plb_acoustics;
// ---------------------------------------------

const T rho0 = 1.;
const T deltaRho = 1.e-4;
const T lattice_speed_sound = 1/sqrt(3);
const T lattice_speed_sound_square = lattice_speed_sound*lattice_speed_sound;
// 120000. 5400 keeps (5000 finish transient and 6000 vortice contact anechoic condition) 
const plint maxIter = 15000;
const plint start_transient_iteration = 5000;
const plint finish_time_to_began_vortice_anechoic_condition = 6000;
//const plint maxIter = 5500; // 120000. 5400 keeps
const plint nx = 1000;       // Choice of lattice dimensions.
const plint ny = 1000;
const T reynolds_number = 150;
const T mach_number = 0.2;
const T velocity_flow = mach_number*lattice_speed_sound;
const plint size_square = 2;
const T tau = (0.5 + ((velocity_flow*size_square)/(reynolds_number*lattice_speed_sound*lattice_speed_sound)));
const T omega = 1/tau;

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;

    pcout << "omega: " << omega << std::endl;
    pcout << "velocity_flow: " << velocity_flow << std::endl;
    pcout << "resultado: " << tau << std::endl;

    pcout << "Total iteration: " << maxIter << std::endl;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new CompleteBGKdynamics<T,DESCRIPTOR>(omega));

    Array<T,2> u0((T)0,(T)0);

    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    
    lattice.initialize();

    plint size_square = 2;
    Box2D square(nx/2 - size_square/2, nx/2 + size_square/2,
    ny/2 - size_square/2, ny/2 + size_square/2);
    defineDynamics(lattice, square, new BounceBack<T,DESCRIPTOR>((T)999));

    // Anechoic Condition
    T rhoBar_target = 0;
    Array<T,2> j_target(velocity_flow, 0.0/std::sqrt(3));
    T size_anechoic_buffer = 30;
    // Define Anechoic Boards
    defineAnechoicBoards(nx, ny, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target,
      rhoBar_target, rhoBar_target, rhoBar_target, rhoBar_target);

    // parameters to FW-HS - Radius 10
    Array<T, 2> center((plint) nx/2, (plint) ny/2);
    plint distance_center = 10;
    FW_H_Surface_square fw_h_surface_square(center, distance_center, maxIter, start_transient_iteration);
    // ------------------

    // Setting pressure points to calculate FFT
    Row pressure_points_partial;
    Row pressure_points_complete;
    // Main loop over time iterations.
    for (plint iT = 0; iT <= maxIter; iT++){

        if (iT >= start_transient_iteration){

            pressure_points_complete.
                push_back((lattice.get(nx/2, (ny/2) + 75*size_square).computeDensity() - rho0)*lattice_speed_sound_square);
            if (iT <= finish_time_to_began_vortice_anechoic_condition){
                pressure_points_partial.
                    push_back((lattice.get(nx/2, (ny/2) + 75*size_square).computeDensity() - rho0)*lattice_speed_sound_square);
            }
        

            fw_h_surface_square.import_pressures_velocities(lattice, iT);
        }

       if (iT%100==0) { 
            pcout << "iT= " << iT << endl;

            /*if (iT>=0){
                ImageWriter<T> imageWriter("leeloo");
                imageWriter.writeScaledGif(createFileName("velocity", iT, 6),
                                   *computeVelocityComponent(lattice, 0));
                imageWriter.writeGif(createFileName("density", iT, 6), 
                *computeDensity(lattice), (T) rho0 + -0.001, (T) rho0 + 0.001); //(T) rho0 + -0.001, (T) rho0 + 0.001);
            }*/
           
            /*
            plb_ofstream matrix_pressure_file("matrix_pressure.dat");
            if (iT == 30000){
                matrix_pressure_file << setprecision(10) << 
            }*/
            
            /*ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(lattice) );*/
                               
           // imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeVorticity(*computeVelocity(lattice)));
            //imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeDensity(lattice));
           //imageWriter.writeGif(createFileName("u", iT, 6), *computeDensity(lattice), (T) rho0 - deltaRho/1000, (T) rho0 + deltaRho/1000);
             pcout << " energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << " rho ="
                  << getStoredAverageDensity<T>(lattice)
                  << " max_velocity ="
                  << setprecision(10) << (getStoredMaxVelocity<T>(lattice))/lattice_speed_sound
                  << endl;
        }

        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
        
    }

    fw_h_surface_square.save_data("sfwh_pressure.dat", "sfwh_velocity_x.dat", "sfwh_velocity_y.dat");

    plb_ofstream pressure_points_complete_file("pressure_points_complete_file.dat");
    plb_ofstream pressure_points_partial_file("pressure_points_complete_partial.dat");
    for (int i = 0; i < pressure_points_complete.size(); ++i){
        pressure_points_complete_file << setprecision(10) << pressure_points_complete[i] << std::endl;
    }
    for (int i = 0; i < pressure_points_partial.size(); ++i){
        pressure_points_partial_file << setprecision(10) << pressure_points_partial[i] << std::endl;
    }

}
