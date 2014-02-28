/*
 * TrackGenerator.h
 *
 *  Created on: Jan 23, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include "Point.h"
#include "Track.h"
#include "Geometry.h"
#include "Plotter.h"


class TrackGenerator {
private:
    /** Number of azimuthal angles in \f$ [0, \pi] \f$ */
    int _num_azim;			

    /** The user-specified track spacing (cm) */
    double _spacing;		

    /** An array of the number of tracks for each azimuthal angle */
    int* _num_tracks;		

    /** An array of the number of tracks starting on the x-axis for each
     *  azimuthal angle */
    int* _num_x;		

    /** An array of the number of tracks starting on the y-axis for each
     *  azimuthal angle */
    int* _num_y;		
    double* _azim_weights;	
    Track** _tracks;
    Geometry* _geom;
    Plotter* _plotter;

    std::string _geometry_file;
    int* _num_segments;
    int _tot_num_tracks;
    int _tot_num_segments;
    Options* _opts;

    /** Boolean or whether to use track input file (true) or not (false) */
    bool _use_input_file;

    /** Filename for the *.tracks input / output file */
    std::string _tracks_filename;

    /** Boolean whether the tracks have been generated (true) or not (false) */
    bool _contains_tracks;
    void computeEndPoint(Point* start, Point* end,  const double phi,
                         const double width, const double height);
    void dumpTracksToFile();
    bool readTracksFromFile();

public:
    TrackGenerator(Geometry* geom, Plotter* plotter, Options* opts);
    virtual ~TrackGenerator();
    double *getAzimWeights() const;
    int getNumAzim() const;
    int *getNumTracks() const;
    double getSpacing() const;
    Track **getTracks() const;
    void generateTracks();

    void makeReflective();
    void segmentize();
    void plotSpec();
    void printTrackingTimers();
    bool containsTracks();

};

#endif /* TRACKGENERATOR_H_ */
