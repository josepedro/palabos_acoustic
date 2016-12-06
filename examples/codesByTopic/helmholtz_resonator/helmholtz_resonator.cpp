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
using namespace plb_acoustics_2D;
// ---------------------------------------------

const T rho0 = 1.;
const T deltaRho = 1.e-2;
const T lattice_speed_sound = 1/sqrt(3);
const T lattice_speed_sound_square = lattice_speed_sound*lattice_speed_sound;
// 120000. 5400 keeps (5000 finish transient and 6000 vortice contact anechoic condition) 
const plint maxIter = 100000;
//const plint maxIter = 5500; // 120000. 5400 keeps
const plint nx = 500;       // Choice of lattice dimensions.
const plint ny = 500;
const T reynolds_number = 150;
const T mach_number = 0.1;
const T velocity_flow = mach_number*lattice_speed_sound;
//const T tau = (0.5 + ((velocity_flow*size_square)/(reynolds_number*lattice_speed_sound*lattice_speed_sound)));
const T omega = 1.98;

T get_linear_chirp(T frequency_min, T frequency_max, plint maxT_final_source, plint iT, T drho){
    T lattice_speed_sound = 1/sqrt(3);
    T initial_frequency = frequency_min;
    T frequency_max_lattice = frequency_max;
    T variation_frequency = (frequency_max_lattice - initial_frequency)/maxT_final_source;
    T frequency_function = initial_frequency*iT + (variation_frequency*iT*iT)/2;
    T phase = 2*M_PI*frequency_function;
    T chirp_hand = 1. + drho*sin(phase);

    return chirp_hand;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;
    pcout << "omega: " << omega << std::endl;
    pcout << "velocity_flow: " << velocity_flow << std::endl;
    pcout << "Total iteration: " << maxIter << std::endl;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new CompleteBGKdynamics<T,DESCRIPTOR>(omega));

    Array<T,2> u0((T)0,(T)0);

    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    
    lattice.initialize();

    // Anechoic Condition
    T rhoBar_target = 0;
    Array<T,2> j_target(0, 0.0/std::sqrt(3));
    T size_anechoic_buffer = 30;
    // Define Anechoic Boards
    defineAnechoicBoards(nx, ny, lattice, size_anechoic_buffer,
      omega, j_target, j_target, j_target, j_target,
      rhoBar_target, rhoBar_target, rhoBar_target, rhoBar_target);

    Box2D above_bounce_back(0, nx - 1, ny - 31 , ny - 1);
    defineDynamics(lattice, above_bounce_back, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_1(0, nx/2 - 100, 0, ny/2);
    defineDynamics(lattice, square_1, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_2(nx/2 - 100, nx/2 + 100, 0, ny/10);
    defineDynamics(lattice, square_2, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_3(nx/2 + 100, nx, 0, ny/2);
    defineDynamics(lattice, square_3, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_4(nx/2 - 100, nx/2 - 10, ny/2 - 10, ny/2);
    defineDynamics(lattice, square_4, new BounceBack<T,DESCRIPTOR>((T)0));

    Box2D square_5(nx/2 + 10, nx/2 + 100, ny/2 - 10, ny/2);
    defineDynamics(lattice, square_5, new BounceBack<T,DESCRIPTOR>((T)0));

    // Files to capture data in time
    Box2D point_capture(nx/2 + 20, nx/2 + 20, ny/2 - 100, ny/2 - 100);
    //Box2D point_capture(nx/2, nx/2, ny/2 - 5, ny/2 - 5);
    plb_ofstream file_pressures("pressures.dat");
    plb_ofstream file_velocities_x("velocities_x.dat");
    plb_ofstream file_velocities_y("velocities_y.dat");

    // Main loop over time iterations.
    Box2D signal_input_place(31, 32, 0, ny - 1);
    for (plint iT = 0; iT <= maxIter; iT++){


        // excitation signal
        if (iT >= 0)
        {
            //T rho_changing = 1. + drho*sin(2*M_PI*(lattice_speed_sound/20)*iT);
            T rho_changing = get_linear_chirp(0, 0.1, maxIter - 1000, iT - 1000, deltaRho);
            initializeAtEquilibrium(lattice, signal_input_place, rho_changing, u0);
        }
        


        if (iT%20==0) { 
            pcout << "iT= " << iT << endl;

            /*VtkImageOutput2D<T> vtkOut(createFileName("vtk", iT, 6), 1.);
            vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
            vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", 1.);
            vtkOut.writeData<float>(*computeDensity(lattice), "vorticity", 1.);*/

            pcout << " energy ="
            << setprecision(10) << getStoredAverageEnergy<T>(lattice)
            << " rho ="
            << getStoredAverageDensity<T>(lattice)
            << " max_velocity ="
            << setprecision(10) << (getStoredMaxVelocity<T>(lattice))///lattice_speed_sound
            << endl;
        }

        file_pressures << setprecision(10) << (computeAverageDensity(lattice, point_capture) - rho0)*lattice_speed_sound_square << endl;
        std::auto_ptr<MultiScalarField2D<T> > velocity_x(plb::computeVelocityComponent(lattice, point_capture, 0));
        file_velocities_x << setprecision(10) << computeAverage(*velocity_x, point_capture) << endl;
        std::auto_ptr<MultiScalarField2D<T> > velocity_y(plb::computeVelocityComponent(lattice, point_capture, 1));
        file_velocities_y << setprecision(10) << computeAverage(*velocity_y, point_capture) << endl;

        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
        
    }

}
