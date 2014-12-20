/*
 * MeshCell.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESHCELL_H_
#define MESHCELL_H_

#define _USE_MATH_DEFINES
#include <assert.h>
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
    double _volume;
    double _at_width;
    double _at_height;
    double _at_volume;
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
    double _diffusion[NUM_ENERGY_GROUPS];
    double _old_flux[NUM_ENERGY_GROUPS];
    double _new_flux[NUM_ENERGY_GROUPS];
    int _cell_id;

    /* for LOO */
    double _l;
    double _at_l;
    double _old_src[NUM_ENERGY_GROUPS]; /* $\bar{Q}^{(m)}$ */
    double _src[NUM_ENERGY_GROUPS]; /* $\bar{Q}^{(m+1/2)}$ */
    double _sum_quad_flux[NUM_ENERGY_GROUPS];
    double _quad_flux[8*NUM_ENERGY_GROUPS];
    double _quad_src[8*NUM_ENERGY_GROUPS];
    double _quad_xs[8*NUM_ENERGY_GROUPS];
    double _quad_in_flux[8*NUM_ENERGY_GROUPS];
    double _quad_out_flux[8*NUM_ENERGY_GROUPS];
    double _net_current[NUM_ENERGY_GROUPS];
    double _old_net_current[NUM_ENERGY_GROUPS];    
public:
    MeshCell();
    virtual ~MeshCell();
    double getLength(int s);
    double getWidth();
    double getATWidth();
    double getHeight();
    double getATHeight();
    double getVolume();
    double getATVolume();
    double getL();
    double getATL();

    void setWidth(double width);
    void setATWidth(double width);
    void setHeight(double height);
    void setATHeight(double height);
    void setVolume(double volume);
    void setATVolume(double volume);
    void setL(double l);
    void setATL(double l);

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
    void setDiffusion(double diffusivity, int e);
    double* getDiffusion();
    MeshSurface* findSurface(LocalCoords* coord, int i);
    void setBounds(double x, double y);
    double* getBounds();
    int getFSRStart();
    int getFSREnd();
    void setFSRStart(int fsr);
    void setFSREnd(int fsr);
    __inline__ MeshSurface* getMeshSurfaces(int surface_id){
        assert(surface_id >= 0);
        assert(surface_id < 8);
        return &_mesh_surfaces[surface_id];
    }
    void setNewFlux(double flux, int e);
    void updateNewFlux(double ratio, int e);
    double* getNewFlux();
    void setOldFlux(double flux, int e);
    double* getOldFlux();
    double getTemp();
    void setTemp(double temp);
    Material* getMaterial();
    void setMaterial(Material* material);
    /* LOO */
    double* getSrc();
    void setSrc(double src, int e);
    double* getOldSrc();
    void setOldSrc(double src, int e);
    double* getQuadFlux();
    double getQuadInFlux(int e, int index);
    double getQuadOutFlux(int e, int index);
    void setQuadFlux(double quadFlux, int e, int index);
    void setQuadInFlux(double quadFlux, int e, int index);
    void setQuadOutFlux(double quadFlux, int e, int index);
    double* getQuadSrc();
    double getQuadSrc(int e, int index);
    void setQuadSrc(double quadSrc, int e, int index);
    double* getQuadXs();
    double getQuadXs(int e, int index);
    void setQuadXs(double quadXs, int e, int index);
    double* getSumQuadFlux();
    void setSumQuadFlux(double sumQuadFlux, int e);
    double* getNetCurrent();
    void setNetCurrent(double current, int e);
    double* getOldNetCurrent();
    void setOldNetCurrent(double current, int e);
};


#endif /* MESHCELL_H_ */
