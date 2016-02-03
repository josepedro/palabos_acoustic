
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

template<typename T, template<typename U> class Descriptor>
void defineAnechoicWall(plint nx, plint ny, MultiBlockLattice2D<T,Descriptor>& lattice, T size_anechoic_buffer, plint side, T omega){
    for(T delta = 0; delta <= size_anechoic_buffer; delta++){        
        DotList2D points_to_aplly_dynamics;
        for (int i = 0; i <= ny; ++i){
            points_to_aplly_dynamics.addDot(Dot2D(100 + delta, i));
        }
        AnechoicDynamics<T,DESCRIPTOR> *anechoicDynamics = new AnechoicDynamics<T,DESCRIPTOR>(omega);
        anechoicDynamics->setDelta((T) delta);
        
        defineDynamics(lattice, points_to_aplly_dynamics, anechoicDynamics);
    }
}