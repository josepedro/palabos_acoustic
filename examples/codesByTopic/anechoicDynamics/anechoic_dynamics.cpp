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
#define PI 3.14159265

// ---------------------------------------------
// Includes of acoustics resources
#include "acoustics/acoustics2D.h"
using namespace plb_acoustics;
typedef vector<T> Row;
typedef vector< Row > Matrix;
// ---------------------------------------------

T rho0 = 1.;
T deltaRho = 1.e-4;
T lattice_speed_sound = 1/sqrt(3);

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;

    const plint maxIter = 15000; // 120000. 5400 keeps
    const plint nx = 1000;       // Choice of lattice dimensions.
    const plint ny = 1000;
    const T omega = 1.98;        // Choice of the relaxation parameter

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
    Array<T,2> j_target(0.15/std::sqrt(3), 0.0/std::sqrt(3));
    T size_anechoic_buffer = 30;
    // Define Anechoic Boards
    defineAnechoicBoards(nx, ny, lattice, size_anechoic_buffer,
	  omega, j_target, j_target, j_target, j_target,
	  rhoBar_target, rhoBar_target, rhoBar_target, rhoBar_target);

    Box2D cima(0, nx, ny - 1, ny);
    //defineDynamics(lattice, cima, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D baixo(0, nx, 0, 1);
    //defineDynamics(lattice, baixo, new BounceBack<T,DESCRIPTOR>(rho0));

    Box2D wall_top(0, nx-1, ny-1, ny-1);
    defineDynamics(lattice, wall_top, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D wall_bottom(0, nx-1, 0, 0);
    defineDynamics(lattice, wall_bottom, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D wall_left(0, 0, 0, ny-1);
    defineDynamics(lattice, wall_left, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D wall_right(nx-1, nx-1, 0, ny-1);
    /*AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
    new AnechoicDynamics<T,DESCRIPTOR>(omega);
    T delta;
    anechoicDynamics->setDelta(delta);
    anechoicDynamics->setRhoBar_target(rhoBar_target);
    anechoicDynamics->setJ_target(j_target);
    //defineDynamics(lattice, wall_top, anechoicDynamics);
    //defineDynamics(lattice, wall_bottom, anechoicDynamics);
    //defineDynamics(lattice, wall_left, anechoicDynamics);
    //defineDynamics(lattice, wall_right, anechoicDynamics);*/

    // parameters to FW-HS
    Array<T, 2> center((plint) nx/2, (plint) ny/2);
    plint distance_center =  200;
    // total number of points
    plint start_transient_iteration = 5400;
    plint total_points_fwhs = distance_center*2*4;
    Matrix matrix_sfwh_pressure(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
    Matrix matrix_sfwh_velocity_x(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
    Matrix matrix_sfwh_velocity_y(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
    // Main loop over time iterations.
    for (plint iT=0; iT <= maxIter; iT++) {

        if (iT >= start_transient_iteration){
        
            plint point_surface = 0;
            // to face 1 (left)
            for (plint y = center[1] - distance_center; y < center[1] + distance_center; y++){
                plint x = center[0] - distance_center;
                matrix_sfwh_pressure[point_surface][0] = x; 
                matrix_sfwh_pressure[point_surface][1] = y;
                matrix_sfwh_velocity_x[point_surface][0] = x; 
                matrix_sfwh_velocity_x[point_surface][1] = y;
                matrix_sfwh_velocity_y[point_surface][0] = x; 
                matrix_sfwh_velocity_y[point_surface][1] = y;
                matrix_sfwh_pressure[point_surface][iT - start_transient_iteration + 2] = (lattice.get(x, y).computeDensity())/3;
                Array<T, 2> velocities((T) 9999, (T) 9999);
                lattice.get(x, y).computeVelocity(velocities);
                matrix_sfwh_velocity_x[point_surface][iT - start_transient_iteration + 2] = velocities[0];
                matrix_sfwh_velocity_y[point_surface][iT - start_transient_iteration + 2] = velocities[1];
                point_surface++;
            }
            // to face 2 (top)
            for (plint x = center[0] - distance_center; x < center[0] + distance_center; x++){
                plint y = center[1] + distance_center;
                matrix_sfwh_pressure[point_surface][0] = x; 
                matrix_sfwh_pressure[point_surface][1] = y;
                matrix_sfwh_velocity_x[point_surface][0] = x; 
                matrix_sfwh_velocity_x[point_surface][1] = y;
                matrix_sfwh_velocity_y[point_surface][0] = x; 
                matrix_sfwh_velocity_y[point_surface][1] = y;
                matrix_sfwh_pressure[point_surface][iT - start_transient_iteration + 2] = (lattice.get(x, y).computeDensity())/3;
                Array<T, 2> velocities((T) 9999, (T) 9999);
                lattice.get(x, y).computeVelocity(velocities);
                matrix_sfwh_velocity_x[point_surface][iT - start_transient_iteration + 2] = velocities[0];
                matrix_sfwh_velocity_y[point_surface][iT - start_transient_iteration + 2] = velocities[1];
                point_surface++;
            }
            // to face 3 (right)
            for (plint y = center[1] + distance_center; y > center[1] - distance_center; y--){
                plint x = center[1] + distance_center;
                matrix_sfwh_pressure[point_surface][0] = x; 
                matrix_sfwh_pressure[point_surface][1] = y;
                matrix_sfwh_velocity_x[point_surface][0] = x; 
                matrix_sfwh_velocity_x[point_surface][1] = y;
                matrix_sfwh_velocity_y[point_surface][0] = x; 
                matrix_sfwh_velocity_y[point_surface][1] = y;
                matrix_sfwh_pressure[point_surface][iT - start_transient_iteration + 2] = (lattice.get(x, y).computeDensity())/3;
                Array<T, 2> velocities((T) 9999, (T) 9999);
                lattice.get(x, y).computeVelocity(velocities);
                matrix_sfwh_velocity_x[point_surface][iT - start_transient_iteration + 2] = velocities[0];
                matrix_sfwh_velocity_y[point_surface][iT - start_transient_iteration + 2] = velocities[1];
                point_surface++;
            }
            // to face 4 (bottom)
            for (plint x = center[0] + distance_center; x > center[0] - distance_center; x--){
                plint y = center[1] - distance_center;
                matrix_sfwh_pressure[point_surface][0] = x; 
                matrix_sfwh_pressure[point_surface][1] = y;
                matrix_sfwh_velocity_x[point_surface][0] = x; 
                matrix_sfwh_velocity_x[point_surface][1] = y;
                matrix_sfwh_velocity_y[point_surface][0] = x; 
                matrix_sfwh_velocity_y[point_surface][1] = y;
                matrix_sfwh_pressure[point_surface][iT - start_transient_iteration + 2] = (lattice.get(x, y).computeDensity())/3;
                Array<T, 2> velocities((T) 9999, (T) 9999);
                lattice.get(x, y).computeVelocity(velocities);
                matrix_sfwh_velocity_x[point_surface][iT - start_transient_iteration + 2] = velocities[0];
                matrix_sfwh_velocity_y[point_surface][iT - start_transient_iteration + 2] = velocities[1];
                point_surface++;
            }

        }

       if (iT%100==0) {  // Write an image every 40th time step.
            pcout << "iT= " << iT << endl;

            /*if (iT>=0){
                ImageWriter<T> imageWriter("leeloo");
                imageWriter.writeScaledGif(createFileName("velocity", iT, 6),
                                   *computeVelocityComponent(lattice, 0));
                imageWriter.writeGif(createFileName("density", iT, 6), 
                *computeDensity(lattice), (T) rho0 + -0.001, (T) rho0 + 0.001); //(T) rho0 + -0.001, (T) rho0 + 0.001);
            }
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
                  << setprecision(10) << getStoredMaxVelocity<T>(lattice)
                  << endl;
        }

        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
        
    }

    plb_ofstream sfwh_pressure_file("data_sfwh/sfwh_pressure.dat");
    plb_ofstream sfwh_velocity_x_file("data_sfwh/sfwh_velocity_x.dat");
    plb_ofstream sfwh_velocity_y_file("data_sfwh/sfwh_velocity_y.dat");
    for (plint point = 0; point < total_points_fwhs; point++){
        for (plint time_step = 0; time_step < maxIter - start_transient_iteration + 2; time_step++){
            sfwh_pressure_file <<  setprecision(10) << matrix_sfwh_pressure[point][time_step] << " ";
            sfwh_velocity_x_file <<  setprecision(10) << matrix_sfwh_velocity_x[point][time_step] << " ";
            sfwh_velocity_y_file <<  setprecision(10) << matrix_sfwh_velocity_y[point][time_step] << " ";
        }
        sfwh_pressure_file << endl;
        sfwh_velocity_x_file << endl;
        sfwh_velocity_y_file << endl;
    }
    sfwh_pressure_file.close();
    sfwh_velocity_x_file.close();
    sfwh_velocity_y_file.close();
}
