/*
 * MeshCell.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESHCELL_H_
#define MESHCELL_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <list>
#include "log.h"
#include <vector>
#include "LocalCoords.h"
#include "MeshSurface.h"

class MeshCell {
private:
	double _width;
	double _height;
	std::vector<int> _FSRs;
	MeshSurface* _mesh_surfaces;
	Material* _material;
	double* _bounds;
	int _fsr_start;
	int _fsr_end;
	double _chi[NUM_ENERGY_GROUPS];
	double _sigma_a[NUM_ENERGY_GROUPS];
	double _sigma_t[NUM_ENERGY_GROUPS];
	double _sigma_s[NUM_ENERGY_GROUPS*NUM_ENERGY_GROUPS];
	double _nu_sigma_f[NUM_ENERGY_GROUPS];
	double _diffusivity[NUM_ENERGY_GROUPS];
	double _old_flux[NUM_ENERGY_GROUPS];
	double _new_flux[NUM_ENERGY_GROUPS];
	double _volume;
	int _cell_id;
	double _temp;
	/* for LOO */
	double _l;
	double _src[NUM_ENERGY_GROUPS]; /* averaged source */
	double _sum_quad_flux[NUM_ENERGY_GROUPS];
	double _quad_flux[8*NUM_ENERGY_GROUPS];
	double _quad_src[8*NUM_ENERGY_GROUPS];

public:
	MeshCell();
	virtual ~MeshCell();
	double getWidth();
	double getHeight();
	double getL();
	void setWidth(double width);
	void setHeight(double height);
	void setL(double l);
	std::vector<int>* getFSRs();
	void addFSR(int fsr);
	void setChi(double chi, int e);
	double* getChi();
	void setSigmaA(double sigmaA, int e);
	double* getSigmaA();
	void setSigmaT(double sigmaT, int e);
	double* getSigmaT();
	void setSigmaS(double sigmaS, int e, int ep);
	double* getSigmaS();
	void setNuSigmaF(double nuSigmaF, int e);
	double* getNuSigmaF();
	void setDiffusivity(double diffusivity, int e);
	double* getDiffusivity();
	MeshSurface* findSurface(LocalCoords* coord, int i);
	void setBounds(double x, double y);
	double* getBounds();
	int getFSRStart();
	int getFSREnd();
	void setFSRStart(int fsr);
	void setFSREnd(int fsr);
	MeshSurface* getMeshSurfaces(int surface);
	void setNewFlux(double flux, int e);
	double* getNewFlux();
	void setOldFlux(double flux, int e);
	double* getOldFlux();
	void setVolume(double volume);
	double getVolume();
	double getTemp();
	void setTemp(double temp);
	Material* getMaterial();
	void setMaterial(Material* material);
	/* LOO */
	void setSrc(double src, int e);
	double* getSrc();
	void setQuadFlux(double quadFlux, int e, int index);
	double* getQuadFlux();
	void setQuadSrc(double quadSrc, int e, int index);
	double* getQuadSrc();
	void setSumQuadFlux(double sumQuadFlux, int e);
	double* getSumQuadFlux();
};


#endif /* MESHCELL_H_ */
