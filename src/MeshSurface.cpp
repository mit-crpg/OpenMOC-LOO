/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){
	
    try{
        _total_wt = new double[6];
        for (int i = 0; i < 6; i++)
            _total_wt[i] = 0.0;
        
        _d_tilde = new double[NUM_ENERGY_GROUPS];
        _d_hat = new double[NUM_ENERGY_GROUPS];
        _current = new double[NUM_ENERGY_GROUPS];
        _d_dif = new double[NUM_ENERGY_GROUPS];
		
        _quad_current = new double*[NUM_ENERGY_GROUPS];
        _quad_flux = new double*[NUM_ENERGY_GROUPS];
        _old_quad_flux = new double*[NUM_ENERGY_GROUPS];
		
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            _quad_current[e] = new double[2];
            _quad_flux[e] = new double[2];
            _old_quad_flux[e] = new double[2];
        }
    }
    catch (std::exception &e)
    {
        log_printf(ERROR, "Unable to allocate memory for mesh surfafces. "
                   "Backtrace:\n%s", e.what());
    }	
		
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
    {
        _d_tilde[e]  = 0.0;
        _d_hat[e]    = 0.0;
        _current[e]  = 0.0;		
        _d_dif[e] = 0.0;		
        /* Assumes 4 quadrature flux per surface, so 2 on each side */
        for (int ind = 0; ind < 2; ind++)
        {
            _quad_current[e][ind] = 0.0;
            _quad_flux[e][ind] = 1.0;
            _old_quad_flux[e][ind] = 1.0;
        }
    }

    _boundary_type = BOUNDARY_NONE;
}

MeshSurface::~MeshSurface()
{
    delete [] _d_tilde;
    delete [] _d_hat;
    delete [] _current;
    delete [] _d_dif;
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
    {
        delete []_quad_current[e];
        delete []_quad_current[e];
    }
    delete []_quad_current;
    delete []_quad_flux;
}

/* Current */
double MeshSurface::getCurrent(int group){
    return _current[group];
}

void MeshSurface::setCurrent(double current, int group){
    _current[group] = current;
}

void MeshSurface::incrementCurrent(double* current){
    for (int group = 0; group < NUM_ENERGY_GROUPS; group++)
    {
        _current[group] += current[group];
    }
}

void MeshSurface::updateCurrent(double factor, int group){
    _current[group] *= factor;
}

/* QuadCurrent */
double MeshSurface::getQuadCurrent(int group, int index){
    return _quad_current[group][index];
}

void MeshSurface::setQuadCurrent(double quad_current, int group, int index){
    _quad_current[group][index] = quad_current;
}

void MeshSurface::incrementQuadCurrent(double quad_current, int group, 
                                       int index)
{	
    _quad_current[group][index] += quad_current;
}
void MeshSurface::incrementTotalWt(double quad_current, int index)
{	
    _total_wt[index] += quad_current;
}

void MeshSurface::setTotalWt(double wt, int index)
{
    _total_wt[index] = wt;
}

void MeshSurface::updateQuadCurrent(double factor, int group, int index)
{	
    _quad_current[group][index] *= factor;
}


/* QuadFlux */
double MeshSurface::getQuadFlux(int group, int index){
    return _quad_flux[group][index];
}

void MeshSurface::setQuadFlux(double quad_flux, int group, int index){
    _quad_flux[group][index] = quad_flux;
}

void MeshSurface::updateQuadFlux(double ratio, int group, int index){
    _quad_flux[group][index] *= ratio;
}

/* OldQuadFlux */
void MeshSurface::setOldQuadFlux(double quad_flux, int group, int index){
    _old_quad_flux[group][index] = quad_flux;
}

double MeshSurface::getOldQuadFlux(int group, int index){
    return _old_quad_flux[group][index];
}



void MeshSurface::setDHat(double dHat, int e){
    _d_hat[e] = dHat;
}

double* MeshSurface::getDHat(){
    return _d_hat;
}

void MeshSurface::setDTilde(double dTilde, int e){
    _d_tilde[e] = dTilde;
}

double* MeshSurface::getDTilde(){
    return _d_tilde;
}

void MeshSurface::setDDif(double dDif, int e){
    _d_dif[e] = dDif;
}

double* MeshSurface::getDDif(){
    return _d_dif;
}

int MeshSurface::getId(){
    return _id;
}

void MeshSurface::setId(int id){
    _id = id;
}

int MeshSurface::getCellId(){
    return _cell_id;
}

void MeshSurface::setCellId(int id){
    _cell_id = id;
}

void MeshSurface::setBoundary(boundaryType boundary){
    _boundary_type = boundary;
}

boundaryType MeshSurface::getBoundary(){
    return _boundary_type;
}



