/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){
	_current = new double[NUM_ENERGY_GROUPS];
	_flux = new double[NUM_ENERGY_GROUPS];

	for (int i = 0; i < NUM_ENERGY_GROUPS; i++){
		_current[i] = 0;
		_flux[i] = 0;
	}
}

MeshSurface::~MeshSurface(){}

meshSurfaceType MeshSurface::getType(){
	return _type;
}

void MeshSurface::setType(meshSurfaceType type){
	_type = type;
}

void MeshSurface::setCurrent(double current, int group){
	_current[group] = current;
}

double MeshSurface::getCurrent(int group){
	return _current[group];
}

void MeshSurface::incrementCurrent(double current, double phi, int group){
	_current[group] += current * fabs(cos(_normal - phi));
}

void MeshSurface::setNormal(double normal){
	_normal = normal;
}

double MeshSurface::getNormal(){
	return _normal;
}

void MeshSurface::setFlux(double flux, int group){
	_flux[group] = flux;
}

double MeshSurface::getFlux(int group){
	return _flux[group];
}

void MeshSurface::incrementFlux(double flux, int group){
	_flux[group] += flux;
}

void MeshSurface::setDHat(double dHat){
	_d_hat = dHat;
}

double MeshSurface::getDHat(){
	return _d_hat;
}

void MeshSurface::setDTilde(double dTilde){
	_d_tilde = dTilde;
}

double MeshSurface::getDTilde(){
	return _d_tilde;
}


