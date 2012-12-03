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
	}
}

MeshSurface::~MeshSurface(){}


void MeshSurface::setCurrent(double current, int group){
	_current[group] = current;
}

double MeshSurface::getCurrent(int group){
	return _current[group];
}

void MeshSurface::incrementCurrent(double flux, int group){
	_current[group] += flux;
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

void MeshSurface::setSurfaceNum(int surfaceNum){
	_surface_num = surfaceNum;
}

int MeshSurface::getSurfaceNum(){
	return _surface_num;
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


