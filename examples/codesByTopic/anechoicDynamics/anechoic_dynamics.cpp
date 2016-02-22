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
T lattice_speed_sound = 1/sqrt(3);

class polar
{
    public:
        float r,th;
        polar(){}
        polar(int a,int b)
        {
            r=a;
            th=b;
        }
        void show()
        {
            cout<<"In polar form:\nr="<<r<<" and theta="<<th;
            
        }
};
class rectangular
{
    public:
        float x,y;
        rectangular(){}
        rectangular(polar p)
        {
            x=p.r*cos(p.th);
            y=p.r*sin(p.th);
        }
        void show()
        {
            cout<<"\nIn Rectangular form:\nx="<<x<<"and y="<<y;
            
        }
};

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;

    const plint maxIter = 120000; // 120000 Iterate during 1000 steps.
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
    defineDynamics(lattice, square, new BounceBack<T,DESCRIPTOR>(rho0));

    // Anechoic Condition
    T rhoBar_target = 0;
    Array<T,2> j_target(0.13/std::sqrt(3), 0.0/std::sqrt(3));
    T size_anechoic_buffer = 30;
  
    //left
    plint orientation = 3;
    Array<T,2> position_anechoic_wall((T)0,(T)0);
    plint length_anechoic_wall = ny + 1;
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall, length_anechoic_wall,
    rhoBar_target, j_target);

    //right
    orientation = 1;
    Array<T,2> position_anechoic_wall_2((T)nx - size_anechoic_buffer,(T)0);
    length_anechoic_wall = ny + 1;
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_2, length_anechoic_wall,
    rhoBar_target, j_target);

    //top
    orientation = 4;
    Array<T,2> position_anechoic_wall_3((T) 0, (T)ny - size_anechoic_buffer);
    length_anechoic_wall = nx + 1;
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_3, length_anechoic_wall,
    rhoBar_target, u0);

    //bottom
    orientation = 2;
    Array<T,2> position_anechoic_wall_1((T)0,(T)0);
    length_anechoic_wall = nx + 1;
    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer, orientation,
    omega, position_anechoic_wall_1, length_anechoic_wall,
    rhoBar_target, u0);

    Box2D cima(0, nx, ny - 1, ny);
    //defineDynamics(lattice, cima, new BounceBack<T,DESCRIPTOR>(rho0));
    Box2D baixo(0, nx, 0, 1);
    //defineDynamics(lattice, baixo, new BounceBack<T,DESCRIPTOR>(rho0));

    // Main loop over time iterations.
    plb_ofstream pressures_time_file("pressures_time.dat");
    plb_ofstream velocities_file("velocities_50000.dat");
    for (plint iT=0; iT<maxIter; ++iT) {
        Box2D centralSquare (nx/2, nx/2, ny/2, ny/2);

        //T rho_changing = 1. + deltaRho*sin(2*PI*(lattice_speed_sound/200)*iT);
        if (iT != 0){
            //initializeAtEquilibrium (lattice, centralSquare, rho_changing, u0);
        }

       if (iT%1000==0) {  // Write an image every 40th time step.
            pcout << "iT= " << iT << endl;

            if (iT>=0){
                ImageWriter<T> imageWriter("leeloo");
                imageWriter.writeScaledGif(createFileName("velocity", iT, 6),
                                   *computeVelocityNorm(lattice) );
                imageWriter.writeGif(createFileName("density", iT, 6), 
                *computeDensity(lattice), (T) rho0 + -0.01, (T) rho0 + 0.01);

                
                // Capturing pressures over time
                T pressure = lattice_speed_sound*lattice_speed_sound*(lattice.get(nx/2, ny - 40).computeDensity() - rho0);
                pressures_time_file << setprecision(10) << pressure << endl;

                // Capturing pressures over space
                plb_ofstream pressures_space_file("pressures_space.dat");
                for (int i = 0; i < ny; i++){
                    T pressure = lattice_speed_sound*lattice_speed_sound*(lattice.get(nx/2, i).computeDensity() - rho0);
                    pressures_space_file << setprecision(10) << pressure << endl;                
                }

                // Capturing pressure matrix
                stringstream ss;
                ss << "matrix_density/density_matrix_" << iT << ".dat";
                //char *name_file = ss.str();
                plb_ofstream pressure_matrix_file(ss.str().c_str());
                pressure_matrix_file << setprecision(10) << *computeDensity(lattice) << endl;
                pressure_matrix_file.close();

            }

            /*plb_ofstream matrix_pressure_file("matrix_pressure.dat");
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

    
   
}
