/*
 * Track.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Track.h"


/*
 * Default track constructor
 */
Track::Track() { }



/**
 * Track destructor
 */
Track::~Track() {
    clearSegments();
#if USE_OPENMP
    omp_destroy_lock(&_flux_lock);
#endif
}


/**
 * Set the values for the track' start and end point and angle
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param phi the track's azimuthal angle
 */
void Track::setValues(const double start_x, const double start_y,
                      const double end_x, const double end_y, const double phi) {

    _start.setCoords(start_x, start_y);
    _end.setCoords(end_x, end_y);
    _phi = phi;

#if USE_OPENMP
    omp_init_lock(&_flux_lock);
#endif
}


/**
 * Set the track azimuthal weight
 * @param weight the azimuthal weight
 */
void Track::setAzimuthalWeight(const double azim_weight) {
    _azim_weight = azim_weight;
}


/**
 * Sets the weight of this track at one of the quadrature polar angles
 * @param angle polar angle
 * @param polar_weight the weight of that angle
 */
void Track::setPolarWeight(const int angle, double polar_weight) {
    _polar_weights[angle] = polar_weight;
}

void Track::updatePolarFluxes(int pe, double update)
{
    _polar_fluxes[pe] *= update;
    return;
}
void Track::setPolarFluxesByIndex(int pe, double flux)
{
    _polar_fluxes[pe] = flux;
    return;
}

void Track::setNewFluxes(int index, double* polar_fluxes)
{
    for (int i = index; i < index + GRP_TIMES_ANG; i++)
        _new_fluxes[i] = polar_fluxes[i];
    return;
}

void Track::zeroNewFluxes()
{
    for (int i = 0; i < 2 * GRP_TIMES_ANG; i++)
        _new_fluxes[i] = 0.0; 
    return;
}

/**
 * Set this track's polar fluxes for a particular direction (0 or 1)
 * @param direction 0 through 3 describing boundaries: REFL_FALSE, REFL_TRUE, 
 *        VAC_FALSE, VAC_TRUE
 * @param polar_fluxes pointer to an array of fluxes
 */
void Track::setPolarFluxes(reflectType direction, double* polar_fluxes) 
{

    if (direction == REFL_TRUE || direction == REFL_FALSE)
    {
        int start = direction * GRP_TIMES_ANG;
        for (int i = 0; i < GRP_TIMES_ANG; i++)
            _polar_fluxes[start + i] = polar_fluxes[i];
    }
    else if (direction == VAC_TRUE || direction == VAC_FALSE)
    {
        int start = (direction - 2) * GRP_TIMES_ANG;
        for (int i = 0; i < GRP_TIMES_ANG; i++)
            _polar_fluxes[start + i] = 0.0;
    }
    else
    {
        log_printf(ERROR, "Tried to set this track's polar flux in a direction"
                   " which does not exist; try using: reflective, vacuume");
    }

    return;
}

/**
 * Set this track's polar fluxes for a particular direction (0 or 1)
 * @param direction 0 through 3 describing boundaries: REFL_FALSE, REFL_TRUE, 
 *        VAC_FALSE, VAC_TRUE
 * @param start_index 0 for forward, and GRP_TIMES_ANG for backward direction
 * @param polar_fluxes pointer to an array of fluxes
 */
void Track::setPolarFluxes(reflectType direction, int start_index,
                           double* polar_fluxes) 
{
    if (direction == REFL_TRUE || direction == REFL_FALSE)
    {
        int start = direction * GRP_TIMES_ANG;
        for (int i = 0; i < GRP_TIMES_ANG; i++)
            _polar_fluxes[start + i] = polar_fluxes[i+start_index];
    }
    else if (direction == VAC_TRUE || direction == VAC_FALSE)
    {
        int start = (direction - 2) * GRP_TIMES_ANG;
        for (int i = 0; i < GRP_TIMES_ANG; i++)
            _polar_fluxes[start + i] = 0.0;
    }
    else
    {
        log_printf(ERROR, "Tried to set this track's polar flux in a direction"
                   " which does not exist; try using: reflective, vacuume");
    }

    return;
}

void Track::setPolarFluxes(reflectType direction, int start_index,
                           double* polar_fluxes, int energy_index) 
{
    if (direction == REFL_TRUE || direction == REFL_FALSE)
    {
        start_index += energy_index * NUM_POLAR_ANGLES;
        int start = direction * GRP_TIMES_ANG + energy_index * NUM_POLAR_ANGLES;
        for (int i = 0; i < NUM_POLAR_ANGLES; i++)
            _polar_fluxes[start + i] = polar_fluxes[start_index + i];
    }
    else if (direction == VAC_TRUE || direction == VAC_FALSE)
    {
        int start = (direction - 2) * GRP_TIMES_ANG 
            + energy_index * NUM_POLAR_ANGLES;
        for (int i = 0; i < NUM_POLAR_ANGLES; i++)
            _polar_fluxes[start + i] = 0.0;
    }
    else
    {
        log_printf(ERROR, "Tried to set this track's polar flux in a direction"
                   " which does not exist; try using: reflective, vacuume");
    }

    return;
}

void Track::printOutInfo()
{
    log_printf(INFO, "(%.3f %.3f)->(%.3f %.3f), forward %f backward %f", 
               _start.getX(), _start.getY(), _end.getX(), _end.getY(),
               _polar_fluxes[0], 
               _polar_fluxes[1]);
}

void Track::printOutNewFluxes()
{
    log_printf(INFO, "(%.3f %.3f)->(%.3f %.3f) %f,"
               " (%.3f %.3f)->(%.3f %.3f) %f",
               _track_out->_start.getX(), _track_out->_start.getY(), 
               _track_out->_end.getX(), _track_out->_end.getY(),
               _new_fluxes[0], 
               _track_in->_start.getX(), _track_in->_start.getY(), 
               _track_in->_end.getX(), _track_in->_end.getY(),        
               _new_fluxes[1]);
}

/**
 * Set this track's polar fluxes for a particular direction (0 or 1)
 * @param direction incoming/outgoing (0/1) flux for forward/reverse directions
 * @param polar_fluxes pointer to an array of fluxes
 */
void Track::resetPolarFluxes(reflectType direction, int start_index)
{
    if (direction == VAC_TRUE || direction == VAC_FALSE)
    {
        int start = (direction - 2) * GRP_TIMES_ANG;
        for (int i = 0; i < GRP_TIMES_ANG; i++)
            _polar_fluxes[start + i] = 0.0;
    }

    return;
}

/*
 * Set the track azimuthal angle
 * @param phi the azimuthal angle
 */
void Track::setPhi(const double phi) {
    _phi = phi;
}


/**
 * Adds a segment pointer to this Track's list of segments
 * IMPORTANT: assumes that segments are added in order of their starting
 * location from the track's start point
 * @param segment a pointer to the segment
 */
void Track::addSegment(segment* segment) {
    try {
        _segments.push_back(segment);
    }
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to add a segment to track. Backtrace:"
                   "\n%s", e.what());
    }
}


/**
 * Sets whether the incoming flux is at the beginning (false) or
 * end (true) of this Track
 * @param relf_in - beginning (false)/end (true)
 */
void Track::setReflIn(reflectType refl_in) {
    _refl_in = refl_in;
}


/**
 * Sets whether the outgoing flux is at the beginning (false) or
 * end (true) of the outgoing Track
 * @param relf_out - beginning (false)/end (true)
 */
void Track::setReflOut(reflectType refl_out) {
    _refl_out = refl_out;
}


/**
 * Sets the incoming track for boundary conditions
 * @param track_in pointer to the incoming track
 */
void Track::setTrackIn(Track *track_in) {
    _track_in = track_in;
}


/**
 * Sets the outgoing track for boundary conditions
 * @param track_out pointer to the outgoing track
 */
void Track::setTrackOut(Track *track_out) {
    _track_out = track_out;
}


/**
 * Return the track's spacing
 * @return the track spacing
 */
void Track::setSpacing(double spacing){
    _spacing = spacing;
}


/**
 * Return the track's azimuthal angle (with respect to the x-axis)
 * @return the azimuthal angle
 */
double Track::getSpacing(){
    return _spacing;
}


/**
 * Returns the track's end point
 * @return a pointer to the track's end point
 */
Point* Track::getEnd() {
    return &_end;
}


/**
 * Returns the track's start point
 * @return a pointer to the track's start point
 */
Point* Track::getStart() {
    return &_start;
}


/**
 * Return the track's azimuthal angle (with respect to the x-axis)
 * @return the aximuthal angle
 */
double Track::getPhi() const {
    return _phi;
}


/**
 * Return the track's azimuthal weight
 * @param the track's azimuthal weight
 */
double Track::getAzimuthalWeight() const {
    return _azim_weight;
}


/**
 * Return an array pointer to the track's polar weights
 * @return pointer to the tracks' polar weights
 */
double* Track::getPolarWeights() {
    return _polar_weights;
}


/**
 * Return a pointer to this track's polar flux array
 * @return a pointer to the polar flux array
 */
double* Track::getPolarFluxes() {
    return _polar_fluxes;
}

double* Track::getNewFluxes() {
    return _new_fluxes;
}

/**
 * Returns the incoming track
 * @return a pointer to the incoming track
 */
Track *Track::getTrackIn() const {
    return _track_in;
}


/**
 * Returns the outgoing track
 * @return a pointer to the outgoing track
 */
Track *Track::getTrackOut() const {
    return _track_out;
}


/**
 * Returns whether the incoming flux is at the start (false) or end
 *  (true) of this Track
 *  @return start (false) or end (true)
 */
reflectType Track::isReflIn() {
    return _refl_in;
}


/**
 * Returns whether the outgoing flux is at the start (false) or end
 * (true) of the outgoing Track
 * @return start (false) or end (true)
 */
reflectType Track::isReflOut() {
    return _refl_out;
}


/**
* Returns a pointer to a segment with a given index or ends program if
* track does not have the segment
* @param segment index into the track's segments container
 * @return a pointer to the requested segment
 */
segment* Track::getSegment(int segment) {

	/* Checks to see if segments container contains this segment index */
	if (segment < (int)_segments.size())
		return _segments.at(segment);

	/* If track doesn't contain this segment, exits program */
	else
            log_printf(ERROR, "Attempted to retrieve segment s = %d but track only"
                       "has %d segments", segment, (int)_segments.size());
	exit(1);
}



/**
 * Returns a vector of this track's segments
 * @return vector of segment pointer
 */
std::vector<segment*> Track::getSegments() {
    return _segments;
}


/**
 * Return the number of segments along this track
 * @return the number of segments
 */
int Track::getNumSegments() {
    return _segments.size();
}


/**
 * Normalizes all of the polar flux values by multiplying by a factor
 * @param factor the factor to scale the flux by
 */
void Track::normalizeFluxes(double factor)  {

#if USE_OPENMP
    omp_set_lock(&_flux_lock);
#endif

    /* Loop over all polar fluxes */
    for (int i = 0; i < 2*GRP_TIMES_ANG; i++)
        _polar_fluxes[i] *= factor;

#if USE_OPENMP
    omp_unset_lock(&_flux_lock);
#endif

    return;
}


/**
 * Checks whether a point is contained along this track
 * @param a pointer to the point of interest
 * @return true if the point is on the track, false otherwise
 */
bool Track::contains(Point* point) {

    double m; 		// the slope of the track

    /* If the point is outside of the bounds of the start and end points of the
     * track it does not lie on the track */
    if (!(((point->getX() <= _start.getX()+1.0E-5 &&
            point->getX() >= _end.getX()-1.0E-5)
           || (point->getX() >= _start.getX()-1.0E-5 &&
               point->getX() <= _end.getX()+1.0E-5))
          &&
          ((point->getY() <= _start.getY()+1.0E-5 &&
            point->getY() >= _end.getY()-1.0E-5)
           || (point->getY() >= _start.getY()-1.0E-5 &&
               point->getY() <= _end.getY()+1.0E-5)))) {

        return false;
    }


    /* If the track is vertical */
    if (fabs(_phi - M_PI / 2) < 1E-10) 
    {
        if (fabs(point->getX() - _start.getX()) < 1E-10)
            return true;
        else
            return false;
    }
    /* If the track is not vertical */
    else 
    {
        m = sin(_phi) / cos(_phi);

        /* Use point-slope formula */
        if (fabs(point->getY() - 
                 (_start.getY() + m * (point->getX() - _start.getX()))) < 1e-10)
            return true;
        else
            return false;
    }
}


/**
 * Deletes each of this track's segments
 */
void Track::clearSegments() 
{
    for (int i=0; i < (int)_segments.size(); i++)
        delete _segments.at(i);

    _segments.clear();
}


/**
 * Convert this track's attributes to a character array
 * @return a character array of this track's attributes
 */
std::string Track::toString() {
    std::stringstream string;
    string << "Track: start, x = " << _start.getX() << ", y = " <<
        _start.getY() << " end, x = " << _end.getX() << ", y = "
           << _end.getY() << ", phi = " << _phi << " azim_weight = " <<
        _azim_weight;

    return string.str();
}


void Track::setSurfFwd(int surfFwd){
    _surf_fwd = surfFwd;
}


int Track::getSurfFwd(){
    return _surf_fwd;
}


void Track::setSurfBwd(int surfBwd){
    _surf_bwd = surfBwd;
}


int Track::getSurfBwd(){
    return _surf_bwd;
}


/**
 * @brief Initializes a track's unique ID. 
 * @details This is set by the trackgenerator to correspond to the track's 
 *          location in a 2D ragged array of all tracks.
 * @param uid the track's unique ID
 */
void Track::setUid(int uid) {
    _uid = uid;
}

/** 
 * @brief Set the index for the track's azimuthal angle index.
 * @details The azimuthal angle index corresponds to a an array of all
 *          azimuthal angles for \f$ \theta \in [0, \pi] \f$ owned by
 *          the TrackGenerator class.
 * @param index the azimuthal angle index
 */
void Track::setAzimAngleIndex(const int index) {
    _azim_angle_index = index;
}


/**
 * @brief Return the index for the track's azimuthal angle (with respect to the
 *        x-axis).
 * @return th azimuthal angle index
 */
int Track::getAzimAngleIndex() const {
    return _azim_angle_index;
}

int Track::getUid() {
    return _uid;
}
