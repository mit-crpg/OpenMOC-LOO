/*
 * MeshCell.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "MeshCell.h"

MeshCell::MeshCell(){
	std::vector<int> _FSRs;
	_volume = 0;
	_bounds = new double[4];
	_cell_id = 0;
	
	/* allocate memory for mesh surfaces */
	_mesh_surfaces = new MeshSurface[8];
	/* set surface id */
	_mesh_surfaces[0].setId(0);		/* left */
	_mesh_surfaces[1].setId(1);		/* bottom */
	_mesh_surfaces[2].setId(2);		/* right */
	_mesh_surfaces[3].setId(3);		/* top */
	_mesh_surfaces[4].setId(4);		/* bottom left */
	_mesh_surfaces[5].setId(5);		/* bottom right */
	_mesh_surfaces[6].setId(6);		/* top right */
	_mesh_surfaces[7].setId(7);		/* top left */


	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
	{
		_chi[e]         = 0.0;
		_nu_sigma_f[e]  = 0.0;
		_sigma_a[e]     = 0.0;
		_sigma_t[e]     = 0.0;
		_diffusivity[e] = 0.0;
		_old_flux[e]    = 0.0;
		_new_flux[e]    = 0.0;
		for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
		{
			_sigma_s[e*NUM_ENERGY_GROUPS + g] = 0.0;
		}
		_src[e] = 0.0;
		for (int k = 0; k < 8; k++) 
		{
			_quad_flux[e * 8 + k] = 0.0;
			_quad_src[e * 8 + k] = 0.0;
		}
	}
}

MeshCell::~MeshCell(){}


double MeshCell::getWidth(){
	return _width;
}

double MeshCell::getHeight(){
	return _height;
}

double MeshCell::getL(){
	return _l;
}

void MeshCell::setWidth(double width){
	_width = width;
}

void MeshCell::setHeight(double height){
	_height = height;
}

void MeshCell::setL(double l){
	_l = l;
}

int MeshCell::getFSRStart(){
	return _fsr_start;
}

int MeshCell::getFSREnd(){
	return _fsr_end;
}

void MeshCell::setFSRStart(int fsr){
	_fsr_start = fsr;
}

void MeshCell::setFSREnd(int fsr){
	_fsr_end = fsr;
}


std::vector<int>* MeshCell::getFSRs(){
	return &_FSRs;

}

void MeshCell::addFSR(int fsr){
	_FSRs.push_back(fsr);
}

void MeshCell::setBounds(double x, double y){
	_bounds[0] = x;
	_bounds[1] = y;
	_bounds[2] = x + _width;
	_bounds[3] = y + _height;
}

double* MeshCell::getBounds(){
	return _bounds;
}

MeshSurface* MeshCell::findSurface(LocalCoords* coord, int i){
	MeshSurface* meshSurface = NULL;
	double x = coord->getX();
	double y = coord->getY();

	/* find which surface coord is on */
	/* left */
	if (fabs(x - _bounds[0]) < 1e-8){
		if (fabs(y - _bounds[1]) > 1e-8 && fabs(y - _bounds[3]) > 1e-8){
			meshSurface = &_mesh_surfaces[0];
		}
		else if (fabs(y - _bounds[3]) < 1e-8){
			meshSurface = &_mesh_surfaces[7];
		}
		else{
			meshSurface = &_mesh_surfaces[4];
		}
	}
	/* right */
	else if (fabs(x - _bounds[2]) < 1e-8){
		if (fabs(y - _bounds[1]) > 1e-8 && fabs(y - _bounds[3]) > 1e-8){
			meshSurface = &_mesh_surfaces[2];
		}
		else if (fabs(y - _bounds[3]) < 1e-8){
			meshSurface = &_mesh_surfaces[6];
		}
		else{
			meshSurface = &_mesh_surfaces[5];
		}
	}
	/* top */
	else if (fabs(y - _bounds[3]) < 1e-8){
		meshSurface = &_mesh_surfaces[3];
	}
	/* bottom */
	else if (fabs(y - _bounds[1]) < 1e-8){
		meshSurface = &_mesh_surfaces[1];
	}

	return meshSurface;
}

MeshSurface* MeshCell::getMeshSurfaces(int surface){
	return &_mesh_surfaces[surface];
}

void MeshCell::setChi(double chi, int e){
	_chi[e] = chi;
}

double* MeshCell::getChi(){
	return _chi;
}

void MeshCell::setSigmaA(double sigmaA, int e){
	_sigma_a[e] = sigmaA;
}

double* MeshCell::getSigmaA(){
	return _sigma_a;
}

void MeshCell::setSigmaT(double sigmaT, int e){
	_sigma_t[e] = sigmaT;
}

double* MeshCell::getSigmaT(){
	return _sigma_t;
}

double* MeshCell::getSigmaS(){
	return _sigma_s;
}

void MeshCell::setSigmaS(double sigmaS, int e, int ep){
	_sigma_s[e*NUM_ENERGY_GROUPS + ep] = sigmaS;
}

double* MeshCell::getNuSigmaF(){
	return _nu_sigma_f;
}

void MeshCell::setNuSigmaF(double nuSigmaF, int e){
	_nu_sigma_f[e] = nuSigmaF;
}

double* MeshCell::getDiffusivity(){
	return _diffusivity;
}

void MeshCell::setDiffusivity(double diffusivity, int e){
	_diffusivity[e] = diffusivity;
}

double* MeshCell::getOldFlux(){
	return _old_flux;
}

void MeshCell::setOldFlux(double flux, int e){
	_old_flux[e] = flux;
}

double* MeshCell::getNewFlux(){
	return _new_flux;
}

void MeshCell::setNewFlux(double flux, int e){
	_new_flux[e] = flux;
}

double* MeshCell::getSrc(){
	return _src;
}

void MeshCell::setSrc(double flux, int e){
	_src[e] = flux;
}

double MeshCell::getVolume(){
	return _volume;
}

void MeshCell::setVolume(double volume){
	_volume = volume;
}

void MeshCell::setTemp(double temp){
	_temp = temp;
}

double MeshCell::getTemp(){
	return _temp;
}

Material* MeshCell::getMaterial(){
	return _material;
}

void MeshCell::setMaterial(Material* material){
	_material = material;
}

void MeshCell::setQuadFlux(double quadFlux, int e, int index){
	_quad_flux[e*NUM_ENERGY_GROUPS + index] = quadFlux;
}

double* MeshCell::getQuadFlux(){
	return _quad_flux;
}

void MeshCell::setQuadSrc(double quadSrc, int e, int index){
	_quad_src[e*NUM_ENERGY_GROUPS + index] = quadSrc;
}

double* MeshCell::getQuadSrc(){
	return _quad_src;
}
