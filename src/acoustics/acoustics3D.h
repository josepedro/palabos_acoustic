
#include "palabos3D.h"
#ifndef PLB_PRECOMPILED // Unless precompiled version is used,
#include "palabos3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

namespace plb_acoustics_3D{

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
	void defineAnechoicBoards(plint nx, plint ny, plint nz,
	 MultiBlockLattice3D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, T omega, 
	  Array<T,3> j_target_normal_z_positive,
	  Array<T,3> j_target_normal_z_negative,
	  Array<T,3> j_target_normal_y_positive,
	  Array<T,3> j_target_normal_y_negative,
	  Array<T,3> j_target_normal_x_positive,
	  Array<T,3> j_target_normal_x_negative,
	  T rhoBar_target){
	  	/* if size_anechoic_buffer is not equal 30, we have to change 
	  	total_distance in file dynamicsTemplates3D.h too*/
	 	for(T delta = 0; delta <= size_anechoic_buffer; delta++){
			// for in all points-cell lattice
	 		for(plint z = 0; z < nz; z++){
				for(plint y = 0; y < ny; y++){
					for(plint x = 0; x < nx; x++){
						// condition to right (1)
						if(x == (nx-delta) &&
						  y >= (delta-1) && 
						  y < (ny-delta) &&
						  z >= (delta-1) && 
						  z < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_x_positive);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
						// condition to bottom (2)
						else if(x >= (delta-1) &&
						  x < (nx-delta) &&
						  y == delta &&
						  z >= (delta-1) && 
						  z < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_y_negative);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
						// condition to left (3)
						else if(x == delta &&
						  y >= (delta-1) &&
						  y < (ny-delta) &&
						  z >= (delta-1) && 
						  z < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_x_negative);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
						// condition to top (4)
						else if(x >= (delta-1) &&
						  x < (nx-delta) &&
						  y == (ny-delta) &&
						  z >= (delta-1) && 
						  z < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_y_positive);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
						// condition to front (5)
						else if(x >= (delta-1) &&
						  x < (nx-delta) &&
						  z == (nz-delta) &&
						  y >= (delta-1) && 
						  y < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_z_positive);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
						// condition to back (6)
						else if(x >= (delta-1) &&
						  x < (nx-delta) &&
						  z == delta &&
						  y >= (delta-1) && 
						  y < (nz-delta)){
							// set delta here
							AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
							new AnechoicDynamics<T,DESCRIPTOR>(omega);
							T delta_efective = 30 - delta;
							anechoicDynamics->setDelta(delta_efective);
							anechoicDynamics->setRhoBar_target(rhoBar_target);
							anechoicDynamics->setJ_target(j_target_normal_z_negative);
							DotList3D points_to_aplly_dynamics;
							points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
							defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
						}
					}
				}
			}
		}
	}

	template<typename T, template<typename U> class Descriptor>
	void defineAnechoicMRTBoards(plint nx, plint ny, plint nz,
	 MultiBlockLattice3D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, T omega, 
	  Array<T,3> j_target_normal_z_positive,
	  Array<T,3> j_target_normal_z_negative,
	  Array<T,3> j_target_normal_y_positive,
	  Array<T,3> j_target_normal_y_negative,
	  Array<T,3> j_target_normal_x_positive,
	  Array<T,3> j_target_normal_x_negative,
	  T rhoBar_target){

		j_target_normal_z_positive = -j_target_normal_z_positive;
	  	j_target_normal_z_negative = -j_target_normal_z_negative;
	  	j_target_normal_y_positive = -j_target_normal_y_positive;
	  	j_target_normal_y_negative = -j_target_normal_y_negative;
	  	j_target_normal_x_positive = -j_target_normal_x_positive;
	  	j_target_normal_x_negative = -j_target_normal_x_negative;

	  	typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;
	  	/* if size_anechoic_buffer is not equal 30, we have to change 
	  	total_distance in file dynamicsTemplates3D.h too*/
	 	for(T delta = 0; delta <= size_anechoic_buffer; delta++){
			// for in all points-cell lattice
	 		for(plint z = 0; z < nz; z++){
				for(plint y = 0; y < ny; y++){
					for(plint x = 0; x < nx; x++){
						// condition to right (1)
                        if(x == (nx-delta) &&
                          y >= (delta-1) && 
                          y < (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_x_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to bottom (2)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          y == delta &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_y_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to left (3)
                        else if(x == delta &&
                          y >= (delta-1) &&
                          y < (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_x_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to top (4)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          y == (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_y_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to front (5)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          z == (nz-delta) &&
                          y >= (delta-1) && 
                          y < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_z_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to back (6)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          z == delta &&
                          y >= (delta-1) && 
                          y < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_z_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
					}
				}
			}
		}
	}

template<typename T, template<typename U> class Descriptor>
	void defineAnechoicMRTBoards_limited(plint nx, plint ny, plint nz,
	 MultiBlockLattice3D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, T omega, 
	  Array<T,3> j_target_normal_z_positive,
	  Array<T,3> j_target_normal_z_negative,
	  Array<T,3> j_target_normal_y_positive,
	  Array<T,3> j_target_normal_y_negative,
	  Array<T,3> j_target_normal_x_positive,
	  Array<T,3> j_target_normal_x_negative,
	  T rhoBar_target, plint off_set_z){

		j_target_normal_z_positive = -j_target_normal_z_positive;
	  	j_target_normal_z_negative = -j_target_normal_z_negative;
	  	j_target_normal_y_positive = -j_target_normal_y_positive;
	  	j_target_normal_y_negative = -j_target_normal_y_negative;
	  	j_target_normal_x_positive = -j_target_normal_x_positive;
	  	j_target_normal_x_negative = -j_target_normal_x_negative;

	  	typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;
	  	/* if size_anechoic_buffer is not equal 30, we have to change 
	  	total_distance in file dynamicsTemplates3D.h too*/
	 	for(T delta = 0; delta <= size_anechoic_buffer; delta++){
			// for in all points-cell lattice
	 		for(plint z = off_set_z; z < nz; z++){
				for(plint y = 0; y < ny; y++){
					for(plint x = 0; x < nx; x++){
						// condition to right (1)
                        if(x == (nx-delta) &&
                          y >= (delta-1) && 
                          y < (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_x_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to bottom (2)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          y == delta &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_y_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to left (3)
                        else if(x == delta &&
                          y >= (delta-1) &&
                          y < (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_x_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to top (4)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          y == (ny-delta) &&
                          z >= (delta-1) && 
                          z < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_y_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to front (5)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          z == (nz-delta) &&
                          y >= (delta-1) && 
                          y < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_z_positive);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
                        // condition to back (6)
                        else if(x >= (delta-1) &&
                          x < (nx-delta) &&
                          z == delta + off_set_z &&
                          y >= (delta-1) && 
                          y < (nz-delta)){
                            // set delta here
                            AnechoicBackgroundDynamics *anechoicDynamics = 
                            new AnechoicBackgroundDynamics(omega);
                            T delta_efective = 30 - delta;
                            anechoicDynamics->setDelta(delta_efective);
                            anechoicDynamics->setRhoBar_target(rhoBar_target);
                            anechoicDynamics->setJ_target(j_target_normal_z_negative);
                            DotList3D points_to_aplly_dynamics;
                            points_to_aplly_dynamics.addDot(Dot3D(x,y,z));
                            defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
                        }
					}
				}
			}
		}
	}


template<typename T, template<typename U> class Descriptor>
	void defineAnechoicMRTWall(plint nx, plint ny, plint nz,
	 MultiBlockLattice3D<T,Descriptor>& lattice,
	  T size_anechoic_buffer, plint orientation, T omega, 
	  Array<plint, 3> position_anechoic_wall, plint length_anechoic_wall,plint width_anechoic_wall,
	  T rhoBar_target, Array<T,3> j_target){

	j_target = -j_target;	
	typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;

		// delta increase to the right
		if(orientation == 1){
			for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		        DotList3D points_to_aplly_dynamics;			// altura em y 

		        for (plint j = 0 ; j <= width_anechoic_wall; ++j){

			        for (int i = 0; i <= length_anechoic_wall; ++i){
			            points_to_aplly_dynamics.addDot(
			            	Dot3D(position_anechoic_wall[0] + delta,
			            	 position_anechoic_wall[1] + i,
			            	 position_anechoic_wall[2]+ j) );
			        }
			    }
		        AnechoicBackgroundDynamics *anechoicDynamics = 
		        new AnechoicBackgroundDynamics(omega);
		        T delta_efective = 30 - delta;
		        anechoicDynamics->setDelta((T) delta_efective);
		        anechoicDynamics->setRhoBar_target(rhoBar_target);
		        //j_target[0] = -j_target[0];  
		        anechoicDynamics->setJ_target(j_target);
		        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
			    
	    	}
		}

		// // delta increase to the left
		// else if(orientation == 3){
		// 	for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
		//         DotList2D points_to_aplly_dynamics;
		//         for (int i = 0; i <= length_anechoic_wall; ++i){
		//             points_to_aplly_dynamics.addDot(
		//             	Dot2D(position_anechoic_wall[0] + delta,
		//             	position_anechoic_wall[1] + i));
		//         }
		//         AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = 
		//         new AnechoicDynamics<T,DESCRIPTOR>(omega);
		//         T delta_left = 30 - delta;
		//         anechoicDynamics->setDelta(delta_left);
		//         anechoicDynamics->setRhoBar_target(rhoBar_target);
		//         anechoicDynamics->setJ_target(j_target);
		//         defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
	 //    	}
		// }

		else{
			cout << "Anechoic Dynamics not Defined." << endl;
			cout << "Choose the correct orientation number." << endl;
		}

	}


}
