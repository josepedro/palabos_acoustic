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

T rho0 = 1.;
T deltaRho = 1.e-4;

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    const plint maxIter = 10*150*sqrt(3); // Iterate during 1000 steps.
    const plint nx = 300;       // Choice of lattice dimensions.
    const plint ny = 300;
    const T omega = 1.98;        // Choice of the relaxation parameter

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );

    lattice.periodicity().toggleAll(false); // Use periodic boundaries.

    Array<T,2> u0((T)0,(T)0);

    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    
    lattice.initialize();

    // Square x1,y1=30,30 and x2,y2=270,270
    Box2D wall_1(30, 30, 30, 270);
    Box2D wall_2(30, 270, 270, 270);
    Box2D wall_3(270, 270, 30, 270);
    Box2D wall_4(30, 270, 30, 30);
    defineDynamics(lattice, wall_1, new BounceBack<T,DESCRIPTOR>);
    defineDynamics(lattice, wall_2, new BounceBack<T,DESCRIPTOR>);
    defineDynamics(lattice, wall_3, new BounceBack<T,DESCRIPTOR>);
    defineDynamics(lattice, wall_4, new BounceBack<T,DESCRIPTOR>);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {

        Box2D centralSquare (150, 150, 150, 150);

        T lattice_speed_sound = 1/sqrt(3);
        T rho_changing = 1. + deltaRho*sin(2*PI*(lattice_speed_sound/20)*iT);
        
        initializeAtEquilibrium (lattice, centralSquare, rho_changing, u0);

        if (iT%40==0) {  // Write an image every 40th time step.
            pcout << "Writing GIF file at iT=" << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            ImageWriter<T> imageWriter("leeloo");
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix

            imageWriter.writeGif(createFileName("u", iT, 6), *computeDensity(lattice), (T) rho0 - deltaRho/1000, (T) rho0 + deltaRho/1000);

            //imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeDensity(lattice));
            cout << setprecision(10) << lattice.get(290, 290).computeDensity() << endl;
        }

        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }

    plint x = 150;
    plb_ofstream ofile("profile.dat");
    for (; x < 300; ++x){
        ofile << setprecision(10) << lattice.get(250, 250).computeDensity() - rho0 << endl;
    }
   
    //Box2D line(150, 300, ny/2, ny/2);
    //ofile << setprecision(3) << *computeDensity(*extractSubDomain(lattice, line)) << endl;
}
