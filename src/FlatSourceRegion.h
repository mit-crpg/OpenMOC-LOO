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
#include "assert.h"

#if USE_OPENMP
#include <omp.h>
#endif

class FlatSourceRegion {
private:
    int _quad_id;
    int _id;
    int _mesh_cell_id;
    Material* _material;
    double _volume;
    double _flux[NUM_ENERGY_GROUPS];
    double _source[NUM_ENERGY_GROUPS];
    double _ratios[NUM_ENERGY_GROUPS];     /* Pre-compute source-over-sigma_t */
#if USE_OPENMP
    omp_lock_t _flux_lock;
#endif
public:
    FlatSourceRegion();
    virtual ~FlatSourceRegion();
    int getQuadId() const;
    int getId() const;
    int getMeshCellId() const;
    Material* getMaterial();
    double getVolume() const;
    double* getFlux();
    double* getSource();
    double* getRatios();
    void setId(int id);
    void setQuadId(int id);
    void setMeshCellId(int id);
    void setMaterial(Material* material);
    void setVolume(double volume);
    void incrementVolume(double volume);
    void setFlux(int energy, double flux);
    void incrementFlux(int energy, double flux);
    void incrementFlux(double* flux);
    void setSource(int energy, double source);
    void normalizeFluxes(double factor);
    void computeRatios();
    void computeRatios(int e);
    double computeFissionRate();
};


#endif /* FLATSOURCEREGION_H_ */
