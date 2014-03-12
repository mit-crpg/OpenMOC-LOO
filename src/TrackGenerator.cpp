#include "TrackGenerator.h"


/**
 * TrackGenerator constructor
 * @param geom a pointer to a geometry object
 * @param plotter a pointer to a plotting object
 * @param num_azim number of azimuthal angles
 * @param spacing track spacing
 */
TrackGenerator::TrackGenerator(Geometry* geom, Plotter* plotter,
                               Options* opts) {
    _plotter = plotter;
    _geom = geom;
    _opts = opts;
    _num_azim = opts->getNumAzim() / 2.0;
    _spacing = opts->getTrackSpacing();
    _geometry_file = opts->getGeometryFile();
    _tot_num_tracks = 0;
    _tot_num_segments = 0;
    _num_segments = NULL;
    _contains_tracks = false;
    _use_input_file = false;
    _tracks_filename = "";
}


/**
 * TrackGenerator destructor frees memory for all tracks
 */
TrackGenerator::~TrackGenerator() {
    delete [] _num_tracks;
    delete [] _num_x;
    delete [] _num_y;
    delete [] _azim_weights;

    for (int i = 0; i < _num_azim; i++)
        delete [] _tracks[i];

    delete [] _tracks;
}


/**
 * Return the azimuthal weights array
 * @return the array of azimuthal weights
 */
double* TrackGenerator::getAzimWeights() const {
    return _azim_weights;
}


/**
 * Return the number of azimuthal angles
 * @return the number of azimuthal angles
 */
int TrackGenerator::getNumAzim() const {
    return _num_azim;
}


/**
 * Return the number of tracks array indexed by azimuthal angle
 * @return the number of tracks
 */
int *TrackGenerator::getNumTracks() const {
    return _num_tracks;
}


/**
 * Return the track spacing array
 * @return the track spacing array
 */
double TrackGenerator::getSpacing() const {
    return _spacing;
}


/**-na 16 -ts 0.2 -ps
 * Return the 2D jagged array of track pointers
 * @return the 2D jagged array of tracks
 */
Track **TrackGenerator::getTracks() const {
    return _tracks;
}


/**
 * Computes the effective angles and track spacings. Computes the number of
 * tracks for each azimuthal angle, allocates memory for all tracks at each
 * angle and sets each track's starting and ending points, azimuthal weight,
 * and azimuthal angle
 */
void TrackGenerator::generateTracks() {
    /* Deletes tracks arrays if tracks have been generated */
    if (_contains_tracks) {
        delete [] _num_tracks;
        delete [] _num_segments;
        delete [] _num_x;
        delete [] _num_y;
        delete [] _azim_weights;

        for (int i = 0; i < _num_azim; i++)
            delete [] _tracks[i];

        delete [] _tracks;
    }

    /** Create tracks directory if one does not yet exist */
    std::stringstream directory;
    directory << getOutputDirectory() << "/tracks";
    struct stat st;
    if (!stat(directory.str().c_str(), &st) == 0)
        mkdir(directory.str().c_str(), S_IRWXU);    

    struct stat buffer;
    std::stringstream test_filename;

    // Manipulates _gometry_file so that there is no "/" symbols. 
    for (unsigned i = 0; i < _geometry_file.length(); i++)
    {
        switch(_geometry_file[i]) {
        case '/':
            _geometry_file[i] = '_';
        }
    }

    test_filename << directory.str() << "/"
                  << _geometry_file << "_"
                  << _num_azim*2.0 << "_"
                  << _spacing;

    if (_opts->getUseUpScatteringXS() == true)
        test_filename << "_up";
    else
        test_filename << "_no"; 
    
    if (_geom->getCmfd() || _geom->getLoo()){
    	test_filename << "_cmfd_"
		      << _geom->getMesh()->getMeshLevel();
    }
    test_filename  << ".data";

    _tracks_filename = test_filename.str();

    if (!stat(_tracks_filename.c_str(), &buffer)) {
        if (readTracksFromFile()) {
            _use_input_file = true;
            _contains_tracks = true;
	}
    }
    
    /* If not tracks input file exists, generate tracks */
    if (_use_input_file == false) 
    {
        /* Allocate memory for the tracks */
        try {
            _num_tracks = new int[_num_azim];
            _num_x = new int[_num_azim];
            _num_y = new int[_num_azim];
            _azim_weights = new double[_num_azim];
            _tracks = new Track*[_num_azim];
        }
        catch (std::exception &e) {
            log_printf(ERROR, "Unable to allocate memory for TrackGenerator. "
                       "Backtrace:\n%s", e.what());
        }

        /* Check to make sure that height, width of the geometry are nonzero */
        if ((_geom->getHeight() <= 0 - ON_SURFACE_THRESH) || 
            (_geom->getWidth() <= 0 - ON_SURFACE_THRESH))
        {
            log_printf(ERROR, "The total height and width of the geometry "
                       "must be nonzero for track generation. Please specify "
                       "the height and width in the geometry input file.");
        }

        try 
        {
            log_printf(NORMAL, "Computing azimuthal angles and track "
                       "spacings...");

            /* Each element in arrays corresponds to a track angle in phi_eff */
            /* Track spacing along x,y-axes, and perpendicular to each track */
            double* dx_eff = new double[_num_azim];
            double* dy_eff = new double[_num_azim];
            double* d_eff = new double[_num_azim];

            /* Effective azimuthal angles with respect to positive x-axis */
            double* phi_eff = new double[_num_azim];

            double x1, x2;
            double iazim = _num_azim*2.0;
            double width = _geom->getWidth();
            double height = _geom->getHeight();

            /* create BitMap for plotting */
            BitMap<int>* bitMap = new BitMap<int>;
            bitMap->pixel_x = _plotter->getBitLengthX();
            bitMap->pixel_y = _plotter->getBitLengthY();
            initialize(bitMap);
            bitMap->geom_x = width;
            bitMap->geom_y = height;
            bitMap->color_type = BLACKWHITE;

            /* Determine azimuthal angles and track spacing */
            for (int i = 0; i < _num_azim; i++) 
            {
                /* desired angle */
                double phi = 2.0 * M_PI / iazim * (0.5 + i);

                /* num intersections with x,y-axes */
                _num_x[i] = (int) (fabs(width / _spacing * sin(phi))) + 1;
                _num_y[i] = (int) (fabs(height / _spacing * cos(phi))) + 1;
                log_printf(DEBUG, "For azimuthal angle %d, num_x = %d,"
                           " num_y = %d", i, _num_x[i], _num_y[i]);

                /* total num of tracks */
                _num_tracks[i] = _num_x[i] + _num_y[i];

                /* effective/actual angle (not the angle we desire, but close */
                phi_eff[i] = atan((height * _num_x[i]) / (width * _num_y[i]));

                /* fix angles in range(pi/2, pi) */
                if (phi > M_PI / 2.0)
                    phi_eff[i] = M_PI - phi_eff[i];

                /* Effective track spacing (not spacing we desire, but close */
                dx_eff[i] = (width / _num_x[i]);
                dy_eff[i] = (height / _num_y[i]);
                d_eff[i] = (dx_eff[i] * sin(phi_eff[i]));
            }

            /* Compute azimuthal weights */
            for (int i = 0; i < _num_azim; i++) 
            {
                if (i < _num_azim - 1)
                    x1 = 0.5 * (phi_eff[i+1] - phi_eff[i]);
                else
                    x1 = 2 * M_PI / 2.0 - phi_eff[i];
                
                if (i >= 1)
                    x2 = 0.5 * (phi_eff[i] - phi_eff[i-1]);
                else
                    x2 = phi_eff[i];

                /* Multiply weight by 2 because angles are in [0, Pi] */
                _azim_weights[i] = (x1 + x2) / (2 * M_PI) * d_eff[i] * 2;
            }

            log_printf(INFO, "Generating track start and end points...");

            /* Compute track starting and end points */
            for (int i = 0; i < _num_azim; i++) 
            {
                /* Tracks for azimuthal angle i */
                _tracks[i] = new Track[_num_tracks[i]];

                /* Compute start points for tracks starting on x-axis */
                for (int j = 0; j < _num_x[i]; j++)
                {
                    _tracks[i][j].getStart()->setCoords(dx_eff[i] * (0.5 + j), 
                                                        0);
                }

                /* Compute start points for tracks starting on y-axis */
                for (int j = 0; j < _num_y[i]; j++) {
                    
                    /* If track points to the upper right */
                    if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) > 0)
                    {
                        _tracks[i][_num_x[i]+j].getStart()
                            ->setCoords(0, dy_eff[i] * (0.5 + j));
                    }
                    /* If track points to the upper left */
                    else if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) < 0)
                    {
                        _tracks[i][_num_x[i]+j].getStart()
                            ->setCoords(width, dy_eff[i] * (0.5 + j));
                    }
                }

                /* Compute the exit points for each track */
                for (int j = 0; j < _num_tracks[i]; j++) 
                {
                    /* Set the track's end point */
                    Point* start = _tracks[i][j].getStart();
                    Point* end = _tracks[i][j].getEnd();
                    computeEndPoint(start, end, phi_eff[i], width, height);

                    /* Set the track's azimuthal weight and spacing*/
                    _tracks[i][j].setAzimuthalWeight(_azim_weights[i]);
                    _tracks[i][j].setPhi(phi_eff[i]);
                    _tracks[i][j].setSpacing(d_eff[i]);
                    
                }
                log_printf(DEBUG, "azimuthal angle = %f", phi_eff[i]);
            }

            /* Recalibrate track start and end points to global origin */
            int uid = 0;
            for (int i = 0; i < _num_azim; i++) 
            {
                _tot_num_tracks += _num_tracks[i];

                for (int j = 0; j < _num_tracks[i]; j++) {

                    _tracks[i][j].setUid(uid);
                    uid++;

                    double x0 = _tracks[i][j].getStart()->getX();
                    double y0 = _tracks[i][j].getStart()->getY();
                    double x1 = _tracks[i][j].getEnd()->getX();
                    double y1 = _tracks[i][j].getEnd()->getY();
                    double new_x0 = x0 - _geom->getWidth()/2.0;
                    double new_y0 = y0 - _geom->getHeight()/2.0;
                    double new_x1 = x1 - _geom->getWidth()/2.0;
                    double new_y1 = y1 - _geom->getHeight()/2.0;
                    double phi = _tracks[i][j].getPhi();
                    _tracks[i][j].setValues(new_x0, new_y0, new_x1, new_y1,phi);
                    _tracks[i][j].setAzimAngleIndex(i);

                    /* Add line to segments bitmap */
                    if (_plotter->plotSpecs() == true)
                        drawLine(bitMap, new_x0, new_y0, new_x1, new_y1, 1);
                }
            }
            if (_plotter->plotSpecs() == true)
                plot(bitMap, "tracks", _plotter->getExtension());
            deleteBitMap(bitMap);

            delete [] dx_eff;
            delete [] dy_eff;
            delete [] d_eff;
            delete [] phi_eff;

            segmentize();
            _contains_tracks = true;
            dumpTracksToFile();
            _use_input_file = true;
        }
        catch (std::exception &e) 
        {
            log_printf(ERROR, "Unable to allocate memory needed to generate "
                       "tracks. Backtrace:\n%s", e.what());
        }
    }

    return;    
}


/**
 * This helper method for generateTracks finds the end point of a given track
 * with a defined start point and an angle from x-axis. This function does not
 * return a value but instead saves the x/y coordinates of the end point
 * directly within the track's end point
 * @param start pointer to the track start point
 * @param end pointer to where the end point should be stored
 * @param phi the azimuthal angle
 * @param width the width of the geometry
 * @param height the height of the geometry
 */
void TrackGenerator::computeEndPoint(Point* start, Point* end,
                                     const double phi, const double width, 
                                     const double height) 
{
    double m = sin(phi) / cos(phi); 		/* slope */
    double yin = start->getY(); 			/* y-coord */
    double xin = start->getX(); 			/* x-coord */

    try 
    {
        Point *points = new Point[4];

        /* Determine all possible points */
        points[0].setCoords(0, yin - m * xin);
        points[1].setCoords(width, yin + m * (width - xin));
        points[2].setCoords(xin - yin / m, 0);
        points[3].setCoords(xin - (yin - height) / m, height);

        /* For each of the possible intersection points */
        for (int i = 0; i < 4; i++) 
        {
            /* neglect the trivial point (xin, yin) */
            if (points[i].getX() == xin && points[i].getY() == yin) { }
            /* The point to return will be within the bounds of the cell */
            else if (points[i].getX() >= 0 && points[i].getX() <= width
                     && points[i].getY() >= 0 && points[i].getY() <= height) 
            {
                end->setCoords(points[i].getX(), points[i].getY());
            }
        }

        delete[] points;
        return;
    }
    catch (std::exception &e) 
    {
        log_printf(ERROR, "Unable to allocate memory for intersection points "
                   "in computeEndPoint TrackGenerator. Backtrace:\n%s",
                   e.what());
    }
}


/**
 *  Implements reflective boundary conditions by setting the incoming
 *  and outgoing tracks for each track using a special indexing scheme
 *  into the trackgenerator's 2d jagged array of tracks
 */
void TrackGenerator::makeReflective() {
    log_printf(INFO, "Creating boundary conditions...");

    int nxi, nyi, nti; /* nx, ny, nt for a particular angle */
    Track *curr;
    Track *refl;

    /* Loop over only half the angles since we will set the pointers for
     * connecting tracks at the same time
     */
    for (int i = 0; i < floor(_num_azim / 2); i++) {
        nxi = _num_x[i];
        nyi = _num_y[i];
        nti = _num_tracks[i];
        curr = _tracks[i];
        refl = _tracks[_num_azim - i - 1];

        /* Loop over all of the tracks for this angle */
        for (int j = 0; j < nti; j++) {
            /* More tracks starting along x-axis than y-axis */
            if (nxi <= nyi) {
                /* Bottom to right hand side */
                if (j < nxi) {
                    curr[j].setTrackIn(&refl[j]);
                    refl[j].setTrackIn(&curr[j]);
                    curr[j].setSurfBwd(3);
                    refl[j].setSurfBwd(3);

                    if (_geom->getSurface(3)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_FALSE);
                        refl[j].setReflIn(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_FALSE);
                        refl[j].setReflIn(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[2 * nxi - 1 - j]);
                    refl[2 * nxi - 1 - j].setTrackIn(&curr[j]);
                    curr[j].setSurfFwd(2);
                    refl[2 * nxi - 1 - j].setSurfBwd(2);

                    if (_geom->getSurface(2)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_FALSE);
                        refl[2 * nxi - 1 - j].setReflIn(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_FALSE);
                        refl[2 * nxi - 1 - j].setReflIn(VAC_TRUE);
                    }
                }

                /* Left hand side to right hand side */
                else if (j < nyi) {

                    curr[j].setTrackIn(&refl[j - nxi]);
                    refl[j - nxi].setTrackOut(&curr[j]);
                    curr[j].setSurfBwd(1);
                    refl[j - nxi].setSurfFwd(1);

                    if (_geom->getSurface(1)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_TRUE);
                        refl[j - nxi].setReflOut(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_TRUE);
                        refl[j - nxi].setReflOut(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[j + nxi]);
                    refl[j + nxi].setTrackIn(&curr[j]);
                    curr[j].setSurfFwd(2);
                    refl[j + nxi].setSurfBwd(2);

                    if (_geom->getSurface(2)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_FALSE);
                        refl[j + nxi].setReflIn(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_FALSE);
                        refl[j + nxi].setReflIn(VAC_TRUE);
                    }
                }

                /* Left hand side to top (j > ny) */
                else {

                    curr[j].setTrackIn(&refl[j - nxi]);
                    refl[j - nxi].setTrackOut(&curr[j]);
                    curr[j].setSurfBwd(1);
                    refl[j - nxi].setSurfFwd(1);


                    if (_geom->getSurface(1)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_TRUE);
                        refl[j - nxi].setReflOut(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_TRUE);
                        refl[j - nxi].setReflOut(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[2 * nti - nxi - j - 1]);
                    refl[2 * nti - nxi - j - 1].setTrackOut(&curr[j]);
                    curr[j].setSurfFwd(4);
                    refl[2 * nti - nxi - j - 1].setSurfFwd(4);

                    if (_geom->getSurface(4)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_TRUE);
                        refl[2 * nti - nxi - j - 1].setReflOut(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_TRUE);
                        refl[2 * nti - nxi - j - 1].setReflOut(VAC_TRUE);
                    }
                }
            }

            /* More tracks starting on y-axis than on x-axis */
            else {
                /* Bottom to top */
                if (j < nxi - nyi) {
                    curr[j].setTrackIn(&refl[j]);
                    refl[j].setTrackIn(&curr[j]);
                    curr[j].setSurfBwd(3);
                    refl[j].setSurfBwd(3);

                    if (_geom->getSurface(3)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_FALSE);
                        refl[j].setReflIn(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_FALSE);
                        refl[j].setReflIn(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[nti - (nxi - nyi) + j]);
                    refl[nti - (nxi - nyi) + j].setTrackOut(&curr[j]);
                    curr[j].setSurfFwd(4);
                    refl[nti - (nxi - nyi) + j].setSurfFwd(4);

                    if (_geom->getSurface(4)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_TRUE);
                        refl[nti - (nxi - nyi) + j].setReflOut(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_TRUE);
                        refl[nti - (nxi - nyi) + j].setReflOut(VAC_TRUE);
                    }
                }

                /* Bottom to right hand side */
                else if (j < nxi) {
                    curr[j].setTrackIn(&refl[j]);
                    refl[j].setTrackIn(&curr[j]);
                    curr[j].setSurfBwd(3);
                    refl[j].setSurfBwd(3);

                    if (_geom->getSurface(3)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_FALSE);
                        refl[j].setReflIn(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_FALSE);
                        refl[j].setReflIn(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[nxi + (nxi - j) - 1]);
                    refl[nxi + (nxi - j) - 1].setTrackIn(&curr[j]);
                    curr[j].setSurfFwd(2);
                    refl[nxi + (nxi - j) - 1].setSurfBwd(2);

                    if (_geom->getSurface(2)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_FALSE);
                        refl[nxi + (nxi - j) - 1].setReflIn(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_FALSE);
                        refl[nxi + (nxi - j) - 1].setReflIn(VAC_TRUE);
                    }
                }

                /* Left-hand side to top (j > nx) */
                else {
                    curr[j].setTrackIn(&refl[j - nxi]);
                    refl[j - nxi].setTrackOut(&curr[j]);
                    curr[j].setSurfBwd(1);
                    refl[j - nxi].setSurfFwd(1);

                    if (_geom->getSurface(1)->getBoundary() == REFLECTIVE){
                        curr[j].setReflIn(REFL_TRUE);
                        refl[j - nxi].setReflOut(REFL_FALSE);
                    }
                    else{
                        curr[j].setReflIn(VAC_TRUE);
                        refl[j - nxi].setReflOut(VAC_FALSE);
                    }

                    curr[j].setTrackOut(&refl[nyi + (nti - j) - 1]);
                    refl[nyi + (nti - j) - 1].setTrackOut(&curr[j]);
                    curr[j].setSurfFwd(4);
                    refl[nyi + (nti - j) - 1].setSurfFwd(4);

                    if (_geom->getSurface(4)->getBoundary() == REFLECTIVE){
                        curr[j].setReflOut(REFL_TRUE);
                        refl[nyi + (nti - j) - 1].setReflOut(REFL_TRUE);
                    }
                    else{
                        curr[j].setReflOut(VAC_TRUE);
                        refl[nyi + (nti - j) - 1].setReflOut(VAC_TRUE);
                    }
                }
            }
        }
    }

    return;
}


/**
 * Generate segments for each track and plot segments
 * in bitmap array.
 */
void TrackGenerator::segmentize() {
    log_printf(NORMAL, "Segmenting tracks...");
    double phi, sin_phi, cos_phi;
    double x0, y0, x1, y1;
    int num_segments;
    Track* track;

    /* create BitMaps for plotting */
    BitMap<int>* bitMapFSR = new BitMap<int>;
    BitMap<int>* bitMap = new BitMap<int>;
    bitMapFSR->pixel_x = _plotter->getBitLengthX();
    bitMapFSR->pixel_y = _plotter->getBitLengthY();
    bitMap->pixel_x = _plotter->getBitLengthX();
    bitMap->pixel_y = _plotter->getBitLengthY();
    initialize(bitMap);
    initialize(bitMapFSR);
    bitMapFSR->geom_x = _geom->getWidth();
    bitMapFSR->geom_y = _geom->getHeight();
    bitMap->geom_x = _geom->getWidth();
    bitMap->geom_y = _geom->getHeight();
    bitMapFSR->color_type = RANDOM;
    bitMap->color_type = RANDOM;

    if (_num_segments != NULL)
        delete [] _num_segments;

    if (!_use_input_file)
    {
#if USE_OPENMP
#pragma omp parallel for private(i, j, k, phi,                  \
                                 sin_phi, cos_phi, track,	\
                                 x0, y0, x1, y1, num_segments)
#endif
        /* Loop over all tracks */
        for (int i = 0; i < _num_azim; i++) {
            phi = _tracks[i][0].getPhi();
            sin_phi = sin(phi);
            cos_phi = cos(phi);

            for (int j = 0; j < _num_tracks[i]; j++){
                track = &_tracks[i][j];
                
                _geom->segmentize(track);

                /* plot segments */
                if (_plotter->plotSpecs() == true)
                {
                    x0 = track->getStart()->getX();
                    y0 = track->getStart()->getY();
                    num_segments = track->getNumSegments();
                    for (int k=0; k < num_segments; k++)
                    {
                        log_printf(DEBUG, "Segmented segment: %i...", k);
                        x1 = x0 + cos_phi * track->getSegment(k)->_length;
                        y1 = y0 + sin_phi * track->getSegment(k)->_length;
                        drawLine(bitMap, x0, y0, x1, y1, 
                                 track->getSegment(k)->_region_id);
                        x0 = x1;
                        y0 = y1;
                    }
                }
                /* done with plotting */
            }
        }
            
        /* Compute the total number of segments in the simulation */
	_num_segments = new int[_tot_num_tracks];
	_tot_num_segments = 0;
	for (int i=0; i < _num_azim; i++) 
        {
	    for (int j=0; j < _num_tracks[i]; j++) 
            {
	        track = &_tracks[i][j];
		_num_segments[track->getUid()] = track->getNumSegments();
	        _tot_num_segments += _num_segments[track->getUid()];
	    }
	}

        log_printf(NORMAL, "Generated %d tracks, %d segments...", 
                   _tot_num_tracks, _tot_num_segments);

        if (_plotter->plotSpecs() == true){
            /* plot segments, FSRs, cells, and materials */
            plot(bitMap, "segments", _plotter->getExtension());
            _plotter->copyFSRMap(bitMapFSR->pixels);
            plot(bitMapFSR, "FSRs", _plotter->getExtension());
            _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                    _geom->getFSRtoCellMap());
            plot(bitMap, "cells", _plotter->getExtension());
            _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                    _geom->getFSRtoMaterialMap());
            plot(bitMap, "materials", _plotter->getExtension());

        }
        /* delete bitMaps */
        deleteBitMap(bitMapFSR);
        deleteBitMap(bitMap);
    }

    return;
}

/**
 * Generate segments for each track and plot segments
 * in bitmap array.
 */
void TrackGenerator::plotSpec() {
    /* create BitMaps for plotting */
    BitMap<int>* bitMapFSR = new BitMap<int>;
    BitMap<int>* bitMap = new BitMap<int>;
    bitMapFSR->pixel_x = _plotter->getBitLengthX();
    bitMapFSR->pixel_y = _plotter->getBitLengthY();
    bitMap->pixel_x = _plotter->getBitLengthX();
    bitMap->pixel_y = _plotter->getBitLengthY();
    initialize(bitMap);
    initialize(bitMapFSR);
    bitMapFSR->geom_x = _geom->getWidth();
    bitMapFSR->geom_y = _geom->getHeight();
    bitMap->geom_x = _geom->getWidth();
    bitMap->geom_y = _geom->getHeight();
    bitMapFSR->color_type = RANDOM;
    bitMap->color_type = RANDOM;

    if (_plotter->plotSpecs() == true){
        _plotter->copyFSRMap(bitMapFSR->pixels);
        plot(bitMapFSR, "FSRs", _plotter->getExtension());
        _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                _geom->getFSRtoCellMap());
        plot(bitMap, "cells", _plotter->getExtension());
        _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                _geom->getFSRtoMaterialMap());
        plot(bitMap, "materials", _plotter->getExtension());
        
    }
    /* delete bitMaps */
    deleteBitMap(bitMapFSR);
    deleteBitMap(bitMap);

    return;
}


/**
 * @brief Writes all track and segment data to a *.tracks file.
 * @details Storing tracks in a file saves time by eliminating segmentation
 *          for commonly simulated geometries.
 */
void TrackGenerator::dumpTracksToFile() {

    if (!_contains_tracks) 
        log_printf(ERROR, "Unable to dump tracks to a file since no tracks "
                 "have been generated for %d azimuthal angles and %f "
                 "track spacing", _num_azim, _spacing);

    FILE* out;
    out = fopen(_tracks_filename.c_str(), "w");

    std::string geometry_to_string = _geom->toString();
    int string_length = geometry_to_string.length();

    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(geometry_to_string.c_str(), sizeof(char)*string_length, 1, out);

    fwrite(&_num_azim, sizeof(int), 1, out);
    fwrite(&_spacing, sizeof(double), 1, out);
    fwrite(_num_tracks, sizeof(int), _num_azim, out);
    fwrite(_num_x, sizeof(int), _num_azim, out);
    fwrite(_num_y, sizeof(int), _num_azim, out);

    double* azim_weights = new double[_num_azim];
    for (int i=0; i < _num_azim; i++)
        azim_weights[i] = _azim_weights[i];
    fwrite(azim_weights, sizeof(double), _num_azim, out);
    delete[] azim_weights;
    
    Track* curr_track;
    double x0, y0, x1, y1;
    double phi;
    int azim_angle_index;
    int num_segments;
    std::vector<segment*> _segments;

    segment* curr_segment;
    double length;
    int material_id;
    int region_id;
    int mesh_surface_fwd;
    int mesh_surface_bwd;

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
            curr_track = &_tracks[i][j];
            x0 = curr_track->getStart()->getX();
            y0 = curr_track->getStart()->getY();
            x1 = curr_track->getEnd()->getX();
            y1 = curr_track->getEnd()->getY();
            phi = curr_track->getPhi();
            azim_angle_index = curr_track->getAzimAngleIndex();
            num_segments = curr_track->getNumSegments();

            fwrite(&x0, sizeof(double), 1, out);
            fwrite(&y0, sizeof(double), 1, out);
            fwrite(&x1, sizeof(double), 1, out);
            fwrite(&y1, sizeof(double), 1, out);
            fwrite(&phi, sizeof(double), 1, out);
            fwrite(&azim_angle_index, sizeof(int), 1, out);
            fwrite(&num_segments, sizeof(int), 1, out);
            
            for (int s=0; s < num_segments; s++) {
                curr_segment = curr_track->getSegment(s);
                length = curr_segment->_length;
                material_id = curr_segment->_material->getId();
                region_id = curr_segment->_region_id;

                fwrite(&length, sizeof(double), 1, out);
                fwrite(&material_id, sizeof(int), 1, out);
                fwrite(&region_id, sizeof(int), 1, out);

                if (_geom->getCmfd() || _geom->getLoo()){
		  mesh_surface_fwd = curr_segment->_mesh_surface_fwd;
		  mesh_surface_bwd = curr_segment->_mesh_surface_bwd;
		  fwrite(&mesh_surface_fwd, sizeof(int), 1, out);
		  fwrite(&mesh_surface_bwd, sizeof(int), 1, out);
		}
            }
        }
    }

    fclose(out);
}


/**
 * @brief Reads tracks in from a file.
 * @details Storing tracks in a file saves time by eliminating segmentation
 *          for commonly simulated geometries.
 */
bool TrackGenerator::readTracksFromFile() {

    /* Deletes tracks arrays if tracks have been generated */
    if (_contains_tracks) 
    {
        delete [] _num_tracks;
        delete [] _num_segments;
        delete [] _num_x;
        delete [] _num_y;
        delete [] _azim_weights;

        for (int i = 0; i < _num_azim; i++)
            delete [] _tracks[i];

        delete [] _tracks;
    }

    FILE* in;
    in = fopen(_tracks_filename.c_str(), "r");

    int string_length;
    int ret = 0;
    ret = fread(&string_length, sizeof(int), 1, in);

    char* geometry_to_string;
    geometry_to_string = new char[string_length+1];

    fread(geometry_to_string, sizeof(char)*string_length, 1, in);
    geometry_to_string[string_length] = '\0';

    if (_geom->toString().compare(0, string_length-1, 
                                  geometry_to_string, 
                                  0, 
                                  string_length-1) != 0) {
        log_printf(NORMAL, "Returning from readTracksFromFile, did not "
                   "find a matching file containing generated tracks.");
        return false;
    }

    log_printf(NORMAL, "Reading segmentized tracks from file %s", 
               _tracks_filename.c_str());

    fread(&_num_azim, sizeof(int), 1, in);
    fread(&_spacing, sizeof(double), 1, in);

    _num_tracks = new int[_num_azim];
    _num_x = new int[_num_azim];
    _num_y = new int[_num_azim];
    _azim_weights = new double[_num_azim];
    double* azim_weights = new double[_num_azim];
    _tracks = new Track*[_num_azim];

    fread(_num_tracks, sizeof(int), _num_azim, in);
    fread(_num_x, sizeof(int), _num_azim, in);
    fread(_num_y, sizeof(int), _num_azim, in);
    fread(azim_weights, sizeof(double), _num_azim, in);

    for (int i = 0; i < _num_azim; i++)
        _azim_weights[i] = azim_weights[i];

    delete[] azim_weights;

    Track* curr_track;
    double x0, y0, x1, y1;
    double phi;
    int azim_angle_index;
    int num_segments;

    double length;
    int material_id;
    int region_id;

    int mesh_surface_fwd;
    int mesh_surface_bwd;
    
    for (int i=0; i < _num_azim; i++)
        _tot_num_tracks += _num_tracks[i];
    _num_segments = new int[_tot_num_tracks];
    int uid = 0;
    _tot_num_segments = 0;
    
    for (int i=0; i < _num_azim; i++) {

        _tracks[i] = new Track[_num_tracks[i]];

        for (int j=0; j < _num_tracks[i]; j++) {

            fread(&x0, sizeof(double), 1, in);
            fread(&y0, sizeof(double), 1, in);
            fread(&x1, sizeof(double), 1, in);
            fread(&y1, sizeof(double), 1, in);
            fread(&phi, sizeof(double), 1, in);
            fread(&azim_angle_index, sizeof(int), 1, in);
            fread(&num_segments, sizeof(int), 1, in);

            _tot_num_segments += num_segments;
            _num_segments[uid] += num_segments;

            curr_track = &_tracks[i][j];
            curr_track->setValues(x0, y0, x1, y1, phi);
            curr_track->setUid(uid);
            curr_track->setAzimAngleIndex(azim_angle_index);
            // Added:
            curr_track->setAzimuthalWeight(_azim_weights[i]);

            for (int s=0; s < num_segments; s++) {
                fread(&length, sizeof(double), 1, in);
                fread(&material_id, sizeof(int), 1, in);
                fread(&region_id, sizeof(int), 1, in);

                segment* curr_segment = new segment;
                curr_segment->_length = length;
                curr_segment->_material = _geom->getMaterial(material_id);
                curr_segment->_region_id = region_id;

                if (_geom->getCmfd() || _geom->getLoo()){
		  fread(&mesh_surface_fwd, sizeof(int), 1, in);
		  fread(&mesh_surface_bwd, sizeof(int), 1, in);
		  curr_segment->_mesh_surface_fwd = mesh_surface_fwd;
		  curr_segment->_mesh_surface_bwd = mesh_surface_bwd;
                }
		
                curr_track->addSegment(curr_segment);
            }

            uid++;
        }
    }

    if (ret)
        _contains_tracks = true;
    fclose(in);    

    delete [] geometry_to_string;
    
    return true;
}
