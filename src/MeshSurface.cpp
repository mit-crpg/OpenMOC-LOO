/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){

	for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
		_d_tilde[e]  = 0.0;
		_d_hat[e]    = 0.0;
		_current[e]  = 0.0;
		
		/* Assumes 4 quadrature flux per surface, so 2 on each side */
		for (int ind = 0; ind < 2; ind++){
			_flux[e][ind] = 0.0;
		}
	}

	_boundary_type = BOUNDARY_NONE;
}

MeshSurface::~MeshSurface(){}

void MeshSurface::setCurrent(double current, int group){
	_current[group] = current;
}

double MeshSurface::getCurrent(int group){
	return _current[group];
}

void MeshSurface::incrementCurrent(double current, int group){
	_current[group] += current;
}

void MeshSurface::setFlux(double flux, int group, int index){
	_flux[group][index] = flux;
}

double MeshSurface::getFlux(int group, int index){
	return _flux[group][index];
}

void MeshSurface::incrementFlux(double flux, int group, int index){
	_flux[group][index] += flux;
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



