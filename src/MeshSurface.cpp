/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){

	_crossings = 0;
	_d_tilde = 0.0;
	_d_hat = 0.0;

}

MeshSurface::~MeshSurface(){}


void MeshSurface::makeCurrents(int numAzim){
	log_printf(DEBUG, "allocating memory for currents and weights");
	_current = new double[numAzim*NUM_ENERGY_GROUPS];
	_current_tot = new double[NUM_ENERGY_GROUPS];
	_weight = new double[numAzim];

}


meshSurfaceType MeshSurface::getType(){
	return _type;
}

void MeshSurface::setType(meshSurfaceType type){
	_type = type;
}

void MeshSurface::setCurrent(double current, int group, int azim){
	_current[azim*NUM_ENERGY_GROUPS + group] = current;
}

double MeshSurface::getCurrent(int group, int azim){
	return _current[azim*NUM_ENERGY_GROUPS + group];
}

void MeshSurface::incrementCurrent(double flux, int group, int azim){
	_current[azim*NUM_ENERGY_GROUPS + group] += flux;
}

void MeshSurface::setCurrentTot(double current, int group){
	_current_tot[group] = current;
}

double MeshSurface::getCurrentTot(int group){
	return _current_tot[group];
}

void MeshSurface::setNormal(double normal){
	_normal = normal;
}

double MeshSurface::getNormal(){
	return _normal;
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

void MeshSurface::setMeshCell(int meshCell){
	_mesh_cell = meshCell;
}

int MeshSurface::getMeshCell(){
	return _mesh_cell;
}

void MeshSurface::setSurfaceNum(int surfaceNum){
	_surface_num = surfaceNum;
}

int MeshSurface::getSurfaceNum(){
	return _surface_num;
}

int MeshSurface::getCrossings(){
	return _crossings;
}

void MeshSurface::incrementCrossings(int crossings){
	_crossings += crossings;
}

void MeshSurface::setCrossings(int crossings){
	_crossings = crossings;
}

double MeshSurface::getWeight(int azim){
	return _weight[azim];
}

void MeshSurface::incrementWeight(double weight, int azim){
	_weight[azim] += weight;
}

void MeshSurface::setWeight(double weight, int azim){
	_weight[azim] = weight;
}

int MeshSurface::getId(){
	return _id;
}

void MeshSurface::setId(int id){
	_id = id;
}






