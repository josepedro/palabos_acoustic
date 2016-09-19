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

/* Orestis Malaspinas contributed this code. */

/** \file
 * MRT dynamics -- generic implementation.
 */
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include "complexDynamics/mrtDynamics.h"
#include "latticeBoltzmann/mrtTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class MRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int MRTdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,MRTdynamics<T,Descriptor> >("MRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>* MRTdynamics<T,Descriptor>::clone() const {
    return new MRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int MRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    //std::cout << "CAI AQUI" << std::endl;
    T jSqr = mrtTemp::mrtCollision(cell, this->getOmega());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class AnechoicMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int AnechoicMRTdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,AnechoicMRTdynamics<T,Descriptor> >("AnechoicMRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AnechoicMRTdynamics<T,Descriptor>::AnechoicMRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AnechoicMRTdynamics<T,Descriptor>* AnechoicMRTdynamics<T,Descriptor>::clone() const {
    return new AnechoicMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int AnechoicMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void AnechoicMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    T jSqr = mrtTemp::anechoicMRTCollision(cell, this->getOmega(),
        this->getDelta(), this->getRhoBar_target(), this->getJ_target());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
void AnechoicMRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T AnechoicMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

// Delta distance because anechoic condition have a buffer
template<typename T, template<typename U> class Descriptor>
void AnechoicMRTdynamics<T,Descriptor>::setDelta(T delta){
    this->delta = delta;
}

template<typename T, template<typename U> class Descriptor>
T AnechoicMRTdynamics<T,Descriptor>::getDelta(){
    return this->delta;
}

// RhoBar_target to develop outflow
template<typename T, template<typename U> class Descriptor>
void AnechoicMRTdynamics<T,Descriptor>::setRhoBar_target(T rhoBar_target){
    this->rhoBar_target = rhoBar_target;
}

template<typename T, template<typename U> class Descriptor>
T AnechoicMRTdynamics<T,Descriptor>::getRhoBar_target(){
    return this->rhoBar_target;
}

// J_target is velocity to anechoic
template<typename T, template<typename U> class Descriptor>
void AnechoicMRTdynamics<T,Descriptor>::setJ_target(Array<T,Descriptor<T>::d> j_target){
    this->j_target = j_target;
}

template<typename T, template<typename U> class Descriptor>
Array<T,Descriptor<T>::d> AnechoicMRTdynamics<T,Descriptor>::getJ_target(){
    return this->j_target;
}

/* *************** Class IncMRTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
int IncMRTdynamics<T,Descriptor>::id =
    meta::registerOneParamDynamics<T,Descriptor,IncMRTdynamics<T,Descriptor> >("IncMRT");
    
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>::IncMRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
IncMRTdynamics<T,Descriptor>* IncMRTdynamics<T,Descriptor>::clone() const {
    return new IncMRTdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
int IncMRTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T jSqr = mrtTemp::incMrtCollision(cell, this->getOmega());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T,Descriptor>::get_rhoBar(cell);
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T thetaBar, BlockStatistics& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;
    
    T jSqr = mrtTemp::incMrtCollision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar) );
    }
}

template<typename T, template<typename U> class Descriptor>
T IncMRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
bool IncMRTdynamics<T,Descriptor>::velIsJ() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void IncMRTdynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell,
                                                                 Array<T,Descriptor<T>::d>& u ) const 
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}


}

#endif  // MRT_DYNAMICS_HH

