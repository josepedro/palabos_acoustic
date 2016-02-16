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
// ---------------------------------------------

T rho0 = 1.;
T deltaRho = 1.e-4;

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;

    const plint maxIter = 200000; // Iterate during 1000 steps.
    const plint nx = 600;       // Choice of lattice dimensions.
    const plint ny = 400;
    const T omega = 1.98;        // Choice of the relaxation parameter

    pcout << "Total iteration: " << maxIter << std::endl;
    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new CompleteBGKdynamics<T,DESCRIPTOR>(omega));

    //lattice.periodicity().toggleAll(false); // Use periodic boundaries.

    Array<T,2> u0((T)0.0/std::sqrt(3),(T)0);

    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    
    lattice.initialize();

    Box2D quadrado( 150 - 10, 150 + 10, 150 - 10, 150 + 10);
    defineDynamics(lattice, quadrado, new BounceBack<T,DESCRIPTOR>(rho0));

    // Anechoic Condition
   T size_anechoic_buffer = 30;
   T rhoBar_target = 1.e-3;
   Array<T,2> j_target(0.11/std::sqrt(3), 0.0/std::sqrt(3));
    //left
   plint orientation = 3;
    Array<T,2> position_anechoic_wall((T)0,(T)0);
    plint length_anechoic_wall = ny + 1;
  defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall, length_anechoic_wall,
    rhoBar_target, j_target);

    //right
    orientation = 1;
    Array<T,2> position_anechoic_wall_2((T)nx - 30,(T)0);
    length_anechoic_wall = ny + 1;
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_2, length_anechoic_wall,
    rhoBar_target, j_target);

    //top
    orientation = 4;
    Array<T,2> position_anechoic_wall_3((T) 0, (T)ny - 30);
    length_anechoic_wall = nx + 1;
    Array<T,2> test(0,0);
    //Array<T,2> test(-j_target[0], j_target[1]);
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_3, length_anechoic_wall,
    rhoBar_target, test);

    //bottom
    orientation = 2;
    Array<T,2> position_anechoic_wall_1((T)0,(T)0);
    length_anechoic_wall = nx + 1;

    //Array<T,2> test_1(-j_target[0], j_target[1]);
    Array<T,2> test_1(0,0);
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_1, length_anechoic_wall,
    rhoBar_target, test_1);

    Box2D cima(0, 300, 298, 299);
    //defineDynamics(lattice, cima, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D baixo(0, 300, 1, 2);
    //defineDynamics(lattice, baixo, new BounceBack<T,DESCRIPTOR>(rho0));


    //defineAnechoicWallOnTheLeftSide(nx, ny, lattice, size_anechoic_buffer, omega);
    //defineAnechoicWallOnTheTopSide(nx, ny, lattice, size_anechoic_buffer, omega);
    //defineAnechoicWallOnTheBottomSide(nx, ny, lattice, size_anechoic_buffer, omega);
    //defineAnechoicWallOnTheRightSide(nx, ny, lattice, size_anechoic_buffer, omega);
    //defineAnechoicWallOnTheLeftSide(nx, ny, lattice, size_anechoic_buffer, omega);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        Box2D centralSquare (150, 150, 150, 150);

        T lattice_speed_sound = 1/sqrt(3);
        T rho_changing = 1. + deltaRho;//*sin(2*PI*(lattice_speed_sound/20)*iT);
        if (iT == 0){
            //initializeAtEquilibrium (lattice, centralSquare, rho_changing, u0);
        }
        
       if (iT%1000==0) {  // Write an image every 40th time step.
            pcout << "iT= " << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            /*ImageWriter<T> imageWriter("leeloo");
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix
            imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(lattice) );*/
           // imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeVorticity(*computeVelocity(lattice)));
            //imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeDensity(lattice));
           //imageWriter.writeGif(createFileName("u", iT, 6), *computeDensity(lattice), (T) rho0 - deltaRho/1000, (T) rho0 + deltaRho/1000);
             pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice)
                  << "; av norm_velocity ="
                  << setprecision(10) << getStoredMaxVelocity<T>(lattice)
                  << endl;
        }

        if (iT == 100){
            plb_ofstream ofile_100("velocity_uy_100.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_100 << setprecision(10) <<  u[1] << endl;    
            }
        }

        if (iT == 500){
            plb_ofstream ofile_500("velocity_uy_500.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_500 << setprecision(10) <<  u[1] << endl;    
            }
        }

        if (iT == 1000){
            plb_ofstream ofile_1000("velocity_uy_1000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_1000 << setprecision(10) <<  u[1] << endl;    
            }
        }

        if (iT == 5000){
            plb_ofstream ofile_5000("velocity_uy_5000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_5000 << setprecision(10) <<  u[1] << endl;    
            }
        }

        if (iT == 10000){
            plb_ofstream ofile_10000("velocity_uy_10000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_10000 << setprecision(10) <<  u[1] << endl;    
            }
        }
        if (iT == 20000){
            plb_ofstream ofile_20000("velocity_uy_20000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_20000 << setprecision(10) <<  u[1] << endl;    
            }
        }
        if (iT == 50000){
            plb_ofstream ofile_50000("velocity_uy_50000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_50000 << setprecision(10) <<  u[1] << endl;    
            }
        }
        if (iT == 100000){
            plb_ofstream ofile_100000("velocity_uy_100000.dat");
            Array<T,2> u;
            for (int i = 0; i < ny; ++i){
                lattice.get(nx/2, i).computeVelocity(u);
                ofile_100000 << setprecision(10) <<  u[1] << endl;    
            }
        }

        
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
        
    }
   
}
