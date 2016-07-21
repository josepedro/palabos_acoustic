
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

namespace plb_acoustics{

	// Defining Matrix and Row
	typedef vector<T> Row;
	typedef vector< Row > Matrix;
	
	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicWall(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, plint orientation, T omega, 
	  Array<plint, 2> position_anechoic_wall, plint length_anechoic_wall,
	  T rhoBar_target, Array<T,2> j_target){
		
		// delta increase to the right
		if(orientation == 1){
			for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		        DotList2D points_to_aplly_dynamics;
		        for (int i = 0; i <= length_anechoic_wall; ++i){
		            points_to_aplly_dynamics.addDot(
		            	Dot2D(position_anechoic_wall[0] + delta,
		            	position_anechoic_wall[1] + i));
		        }
		        AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
		        new AnechoicDynamics<T,DESCRIPTOR>(omega);
		        anechoicDynamics->setDelta((T) delta);
		        anechoicDynamics->setRhoBar_target(rhoBar_target);
		        //j_target[0] = -j_target[0];  
		        anechoicDynamics->setJ_target(j_target);
		        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
	    	}
		}

		// delta increase to the top
		else if(orientation == 4){
			for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		        DotList2D points_to_aplly_dynamics;
		        for (int i = 0; i <= length_anechoic_wall; ++i){
		            points_to_aplly_dynamics.addDot(
		            	Dot2D(position_anechoic_wall[0] + i,
		            	position_anechoic_wall[1] + delta));
		        }
		        AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
		        new AnechoicDynamics<T,DESCRIPTOR>(omega);
		        anechoicDynamics->setDelta((T) delta);
		        anechoicDynamics->setRhoBar_target(rhoBar_target);
		        anechoicDynamics->setJ_target(j_target);
		        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
	    	}
		}

		// delta increase to the left
		else if(orientation == 3){
			for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		        DotList2D points_to_aplly_dynamics;
		        for (int i = 0; i <= length_anechoic_wall; ++i){
		            points_to_aplly_dynamics.addDot(
		            	Dot2D(position_anechoic_wall[0] + delta,
		            	position_anechoic_wall[1] + i));
		        }
		        AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
		        new AnechoicDynamics<T,DESCRIPTOR>(omega);
		        T delta_left = 30 - delta;
		        anechoicDynamics->setDelta(delta_left);
		        anechoicDynamics->setRhoBar_target(rhoBar_target);
		        anechoicDynamics->setJ_target(j_target);
		        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
	    	}
		}

		// delta increase to the bottom
		else if(orientation == 2){
			for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		        DotList2D points_to_aplly_dynamics;
		        for (int i = 0; i <= length_anechoic_wall; ++i){
		            points_to_aplly_dynamics.addDot(
		            	Dot2D(position_anechoic_wall[0] + i,
		            	position_anechoic_wall[1] +  delta));
		        }
		        AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
		        new AnechoicDynamics<T,DESCRIPTOR>(omega);
		        T delta_left = 30 - delta;
		        anechoicDynamics->setDelta(delta_left);
		        anechoicDynamics->setRhoBar_target(rhoBar_target);
		        anechoicDynamics->setJ_target(j_target);
		        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
	    	}
		}

		else{
			cout << "Anechoic Dynamics not Defined." << endl;
			cout << "Choose the correct orientation number." << endl;
		}

	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicBoards(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, T omega, Array<T,2> j_target_1, 
	  Array<T,2> j_target_2, Array<T,2> j_target_3,
          Array<T,2> j_target_4, T rhoBar_target_1, 
	  T rhoBar_target_2, T rhoBar_target_3, T rhoBar_target_4){
	 	for(T delta = 0; delta <= size_anechoic_buffer; delta++){
			// for in all points-cell lattice
			for(plint y = 0; y < ny; y++){
				for(plint x = 0; x < nx; x++){
					// condition to right (1)
					if(x == (nx-delta) &&
					  y >= (delta-1) && 
					  y < (ny-delta)){
						// set delta here
						AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
						new AnechoicDynamics<T,DESCRIPTOR>(omega);
						T delta_efective = 30 - delta;
						anechoicDynamics->setDelta(delta_efective);
						anechoicDynamics->setRhoBar_target(rhoBar_target_1);
						anechoicDynamics->setJ_target(j_target_1);
						DotList2D points_to_aplly_dynamics;
						points_to_aplly_dynamics.addDot(Dot2D(x,y));
						defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
					}
					// condition to bottom (2)
					else if(x >= (delta-1) &&
					  x < (nx-delta) &&
					  y == delta){
						// set delta here
						AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
						new AnechoicDynamics<T,DESCRIPTOR>(omega);
						T delta_efective = 30 - delta;
						anechoicDynamics->setDelta(delta_efective);
						anechoicDynamics->setRhoBar_target(rhoBar_target_2);
						anechoicDynamics->setJ_target(j_target_2);
						DotList2D points_to_aplly_dynamics;
						points_to_aplly_dynamics.addDot(Dot2D(x,y));
						defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
					}
					// condition to left (3)
					else if(x == delta &&
					  y >= (delta-1) &&
					  y < (ny-delta)){
						// set delta here
						AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
						new AnechoicDynamics<T,DESCRIPTOR>(omega);
						T delta_efective = 30 - delta;
						anechoicDynamics->setDelta(delta_efective);
						anechoicDynamics->setRhoBar_target(rhoBar_target_3);
						anechoicDynamics->setJ_target(j_target_3);
						DotList2D points_to_aplly_dynamics;
						points_to_aplly_dynamics.addDot(Dot2D(x,y));
						defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
					}
					// condition to top (4)
					else if(x >= (delta-1) &&
					  x < (nx-delta) &&
					  y == (ny-delta)){
						// set delta here
						AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
						new AnechoicDynamics<T,DESCRIPTOR>(omega);
						T delta_efective = 30 - delta;
						anechoicDynamics->setDelta(delta_efective);
						anechoicDynamics->setRhoBar_target(rhoBar_target_4);
						anechoicDynamics->setJ_target(j_target_4);
						DotList2D points_to_aplly_dynamics;
						points_to_aplly_dynamics.addDot(Dot2D(x,y));
						defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
					}
				}
			}
		}
	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicWallOnTheRightSide(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice, T size_anechoic_buffer, T omega){
	 	plint orientation = 1;
	    plint length_anechoic_wall = ny;
	    // position x and y
	    Array<plint, 2> position_anechoic_wall((plint) nx - size_anechoic_buffer - 1, 0);
	    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer,
	                       orientation, omega, position_anechoic_wall, length_anechoic_wall);
	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicWallOnTheLeftSide(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice, T size_anechoic_buffer, T omega){
	 	plint orientation = 3;
	    plint length_anechoic_wall = ny;
	    // position x and y
	    Array<plint, 2> position_anechoic_wall(0, 0);
	    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer,
	                       orientation, omega, position_anechoic_wall, length_anechoic_wall);
	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicWallOnTheTopSide(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice, T size_anechoic_buffer, T omega){
	 	plint orientation = 4;
	    plint length_anechoic_wall = nx;
	    // position x and y
	    Array<plint, 2> position_anechoic_wall(0, (plint) ny - size_anechoic_buffer - 1);
	    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer,
	                       orientation, omega, position_anechoic_wall, length_anechoic_wall);
	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicWallOnTheBottomSide(plint nx, plint ny,
	 MultiBlockLattice2D<T,Descriptor>& lattice, T size_anechoic_buffer, T omega){
	 	plint orientation = 2;
	    plint length_anechoic_wall = nx;
	    // position x and y
	    Array<plint, 2> position_anechoic_wall(0, 0);
	    defineAnechoicWall(nx, ny, lattice, size_anechoic_buffer,
	                       orientation, omega, position_anechoic_wall, length_anechoic_wall);
	}

	class FW_H_Surface_square {
	    Array<T, 2> center;
	    plint radius;
	    Matrix matrix_sfwh_pressure;
	    Matrix matrix_sfwh_velocity_x;
	    Matrix matrix_sfwh_velocity_y;
	    plint total_points_fwhs;
	    plint maxIter;
	    plint start_transient_iteration;

	  	public:
	  	FW_H_Surface_square(Array<T, 2> center, plint radius, plint maxIter, plint start_transient_iteration);
	  	void import_pressures_velocities(MultiBlockLattice2D<T, DESCRIPTOR> lattice, plint iT);
	  	void save_data(char *pressure_file_name, char *velocity_x_file_name, char *velocity_y_file_name);
	};

	FW_H_Surface_square::FW_H_Surface_square(Array<T,2> center, plint radius, plint maxIter, plint start_transient_iteration){
		this->center = center;
		this->radius = radius;
		this->total_points_fwhs = radius*2*4;
		this->maxIter = maxIter;
		this->start_transient_iteration = start_transient_iteration;
		Matrix matrix_sfwh_pressure(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
		this->matrix_sfwh_pressure = matrix_sfwh_pressure;
    	Matrix matrix_sfwh_velocity_x(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
    	this->matrix_sfwh_velocity_x = matrix_sfwh_velocity_x;
    	Matrix matrix_sfwh_velocity_y(total_points_fwhs, Row(maxIter - start_transient_iteration + 2));
    	this->matrix_sfwh_velocity_y = matrix_sfwh_velocity_y;
		pcout << "raio dessa parada: " << this->radius << std::endl;
	}

	void FW_H_Surface_square::import_pressures_velocities(MultiBlockLattice2D<T, DESCRIPTOR> lattice, plint iT){
		plint point_surface = 0;
        // to face 1 (left)
        for (plint y = this->center[1] - this->radius; y < this->center[1] + this->radius; y++){
            plint x = this->center[0] - this->radius;
            matrix_sfwh_pressure[point_surface][0] = x; 
            matrix_sfwh_pressure[point_surface][1] = y;
            matrix_sfwh_velocity_x[point_surface][0] = x; 
            matrix_sfwh_velocity_x[point_surface][1] = y;
            matrix_sfwh_velocity_y[point_surface][0] = x; 
            matrix_sfwh_velocity_y[point_surface][1] = y;
            matrix_sfwh_pressure[point_surface][iT - start_transient_iteration + 2] = 
                (lattice.get(x, y).computeDensity())/3;
            Array<T, 2> velocities((T) 9999, (T) 9999);
            lattice.get(x, y).computeVelocity(velocities);
            matrix_sfwh_velocity_x[point_surface][iT - start_transient_iteration + 2] = velocities[0];
            matrix_sfwh_velocity_y[point_surface][iT - start_transient_iteration + 2] = velocities[1];
            point_surface++;
        }
        // to face 2 (top)
        for (plint x = this->center[0] - this->radius; x < this->center[0] + this->radius; x++){
            plint y = this->center[1] + this->radius;
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
        for (plint y = this->center[1] + this->radius; y > this->center[1] - this->radius; y--){
            plint x = this->center[1] + this->radius;
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
        for (plint x = this->center[0] + this->radius; x > this->center[0] - this->radius; x--){
            plint y = this->center[1] - this->radius;
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

	void FW_H_Surface_square::save_data(char *pressure_file_name, 
		char *velocity_x_file_name, char *velocity_y_file_name){
		plb_ofstream sfwh_pressure_file(pressure_file_name);
	    plb_ofstream sfwh_velocity_x_file(velocity_x_file_name);
	    plb_ofstream sfwh_velocity_y_file(velocity_y_file_name);
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

}
