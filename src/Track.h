/*
 * Track.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TRACK_H_
#define TRACK_H_

#include <vector>
#include <stdlib.h>
#include <string>
#include "Point.h"
#include "Material.h"
#include "log.h"
#include "configurations.h"

#if USE_OPENMP
#include <omp.h>
#endif

/**
 * Surface boundary types
 */
enum reflectType {
    REFL_FALSE,
    REFL_TRUE,
    VAC_FALSE,
    VAC_TRUE
};



/* Represent a segment along a given track */
struct segment {
    double _length;
    Material* _material;
    int _region_id;
#if STORE_PREFACTORS
    double _prefactors[NUM_ENERGY_GROUPS][NUM_POLAR_ANGLES];
#endif
    int _mesh_surface_fwd;
    int _mesh_surface_bwd;
};

class Track {
private:
    Point _start;
    Point _end;
    double _phi;
    double _spacing;
    double _azim_weight;
    double _polar_weights[NUM_POLAR_ANGLES];
    double _polar_fluxes[2*GRP_TIMES_ANG];
    double _bwd_fluxes[2*GRP_TIMES_ANG];
    double _fwd_fluxes[2*GRP_TIMES_ANG];
    std::vector<segment*> _segments;

    /** The track which reflects out of this track along its "forward"
     * direction for reflective boundary conditions. */
    Track* _track_in;

    /** The track which reflects out of this track along its "reverse"
     * direction for reflective boundary conditions. */
    Track* _track_out;

    /* whether to give the flux to the forward (false) or backward
     * (true) direction of the track reflecting out of this one along
     * its forward (in) or backward (out) direction. */
    reflectType _refl_in, _refl_out;
    int _surf_fwd;
    int _surf_bwd;
    /** A monotonically increasing unique ID for each track */
    int _uid;
    /** The azimuthal angle index into the global 2D ragged array of tracks */
    int _azim_angle_index;
#if USE_OPENMP
    omp_lock_t _flux_lock;
#endif
public:
    Track();
    virtual ~Track();
    void setValues(const double start_x, const double start_y,
                   const double end_x, const double end_y, const double phi);
    void setAzimuthalWeight(const double azim_weight);
    void setPolarWeight(const int angle, double polar_weight);
    void setPolarFluxes(reflectType direction, int start_index, 
                        double* polar_fluxes);
    void resetPolarFluxes(reflectType direction, int start_index);
    void setPolarFluxesByIndex(int pe, double flux);
    void updatePolarFluxes(int pe, double factor);
    void setUid(int uid);
    void setAzimAngleIndex(const int index);
    void setPhi(const double phi);
    void setReflIn(reflectType refl_in);
    void setReflOut(reflectType refl_out);
    void setTrackIn(Track *track_in);
    void setTrackOut(Track *track_out);
    void setSpacing(double spacing);
    double getSpacing();

    Point* getEnd();
    Point* getStart();
    int getUid();
    double getPhi() const;
    int getAzimAngleIndex() const;
    double getAzimuthalWeight() const;
    double* getPolarWeights();
    double* getPolarFluxes();
    double* getFwdFluxes();
    double* getBwdFluxes();
    double* getNewPolarFluxes();
    segment* getSegment(int s);
    std::vector<segment*> getSegments();
    int getNumSegments();
    Track *getTrackIn() const;
    Track *getTrackOut() const;
    reflectType isReflIn();
    reflectType isReflOut();
    void setSurfFwd(int surfFwd);
    int getSurfFwd();
    void setSurfBwd(int surfBwd);
    int getSurfBwd();

    void normalizeFluxes(double factor);
    bool contains(Point* point);
    void addSegment(segment* segment);
    void clearSegments();
    std::string toString();

    void setFwdFluxes(double *fluxes);
    void setBwdFluxes(double *fluxes);

};

#endif /* TRACK_H_ */
