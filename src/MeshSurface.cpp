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
	_crossings = 0;
	_weight = 0;
	_curr_weight = 0;

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

void MeshSurface::incrementCurrent(double flux, double phi, int group){
	_current[group] += flux * fabs(cos(_normal - phi));
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

double MeshSurface::getWeight(){
	return _weight;
}

void MeshSurface::incrementWeight(double weight){
	_weight += weight;
}

void MeshSurface::setWeight(double weight){
	_weight = weight;
}

double MeshSurface::getCurrWeight(){
	return _curr_weight;
}

void MeshSurface::incrementCurrWeight(double weight, double phi){
	_curr_weight += weight * fabs(cos(_normal - phi));
}

void MeshSurface::setCurrWeight(double weight){
	_curr_weight = weight;
}







