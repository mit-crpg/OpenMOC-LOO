/*
 * FlatSourceRegion.h
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef FLATSOURCEREGION_H_
#define FLATSOURCEREGION_H_

#include "Material.h"
#include "configurations.h"
#include "log.h"

#if USE_OPENMP
#include <omp.h>
#endif

class FlatSourceRegion {
private:
    int _id;
    Material* _material;
    double _volume;
    double _flux[NUM_ENERGY_GROUPS];
    double _old_flux[NUM_ENERGY_GROUPS];
    double _source[NUM_ENERGY_GROUPS];
    double _old_source[NUM_ENERGY_GROUPS];
    /* Pre-computed Ratio of source / sigma_t */
    double _ratios[NUM_ENERGY_GROUPS];
    //double _boundary_update[NUM_ENERGY_GROUPS * 2];
#if USE_OPENMP
    omp_lock_t _flux_lock;
#endif
public:
    FlatSourceRegion();
    virtual ~FlatSourceRegion();
    int getId() const;
    Material* getMaterial();
    double getVolume() const;
    double* getFlux();
    double* getOldFlux();
    double* getOldSource();
    double* getSource();
    double* getRatios();
    void setId(int id);
    void setMaterial(Material* material);
    void setVolume(double volume);
    void incrementVolume(double volume);
    void setFlux(int energy, double flux);
    void incrementFlux(int energy, double flux);
    void incrementFlux(double* flux);
    void setOldFlux(int energy, double old_flux);
    void setSource(int energy, double source);
    void setOldSource(int energy, double old_source);
    void normalizeFluxes(double factor);
    void computeRatios();
    double computeFissionRate();
    /*
    void setBoundaryUpdate(int group, int ind, double boundary_update);
    double getBoundaryUpdate(int group, int ind);
    */
};


#endif /* FLATSOURCEREGION_H_ */
