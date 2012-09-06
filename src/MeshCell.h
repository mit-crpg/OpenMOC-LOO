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
	double _abs_rate;
	MeshSurface* _mesh_surfaces;
	double* _bounds;
	int _fsr_start;
	int _fsr_end;
	double _sigma_a;
	double _sigma_t;
	double _sigma_f;
	double _sigma_s;
	double _nu_sigma_f;
	double _diffusivity;
	double _old_flux;
	double _new_flux;

public:
	MeshCell();
	virtual ~MeshCell();
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	std::vector<int>* getFSRs();
	void addFSR(int fsr);
	void setAbsRate(double absRate);
	double getAbsRate();
	void setSigmaA(double sigmaA);
	double getSigmaA();
	void setSigmaT(double sigmaT);
	double getSigmaT();
	void setSigmaF(double sigmaF);
	double getSigmaF();
	void setSigmaS(double sigmaS);
	double getSigmaS();
	void setNuSigmaF(double nuSigmaF);
	double getNuSigmaF();
	void setDiffusivity(double diffusivity);
	double getDiffusivity();
	MeshSurface* findSurface(LocalCoords* coord);
	void setBounds(double x, double y);
	double* getBounds();
	int getFSRStart();
	int getFSREnd();
	void setFSRStart(int fsr);
	void setFSREnd(int fsr);
	void makeSurfaces();
	MeshSurface* getMeshSurfaces(int surface);
	void setNewFlux(double flux);
	double getNewFlux();
	void setOldFlux(double flux);
	double getOldFlux();
};


#endif /* MESHCELL_H_ */
