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


/**
 * Set this track's polar fluxes for a particular direction (0 or 1)
 * @param direction incoming/outgoing (0/1) flux for forward/reverse directions
 * @param polar_fluxes pointer to an array of fluxes
 */
void Track::setPolarFluxes(bool direction, int start_index,
							double* polar_fluxes) {
#if USE_OPENMP
	omp_set_lock(&_flux_lock);
#endif

	int start = direction * GRP_TIMES_ANG;

	if (direction != true && direction != false)
		log_printf(ERROR, "Tried to set this track's polar flux in a direction"
				"which does not exist: direction = %b", direction);

	for (int i = 0; i < GRP_TIMES_ANG; i++)
		_polar_fluxes[start + i] = polar_fluxes[i+start_index];

#if USE_OPENMP
	omp_unset_lock(&_flux_lock);
#endif

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
void Track::setReflIn(const bool refl_in) {
    _refl_in = refl_in;
}


/**
 * Sets whether the outgoing flux is at the beginning (false) or
 * end (true) of the outgoing Track
 * @param relf_out - beginning (false)/end (true)
 */
void Track::setReflOut(const bool refl_out) {
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
bool Track::isReflIn() const {
    return _refl_in;
}


/**
 * Returns whether the outgoing flux is at the start (false) or end
 * (true) of the outgoing Track
 * @return start (false) or end (true)
 */
bool Track::isReflOut() const {
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
				"has %d segments", segment, _segments.size());
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
	if (!(((point->getX() <= _start.getX()+1.0E-2 &&
			point->getX() >= _end.getX()-1.0E-2)
		|| (point->getX() >= _start.getX()-1.0E-2 &&
				point->getX() <= _end.getX()+1.0E-2))
		&&
		((point->getY() <= _start.getY()+1.0E-2 &&
				point->getY() >= _end.getY()-1.0E-2)
		|| (point->getY() >= _start.getY()-1.0E-2 &&
				point->getY() <= _end.getY()+1.0E-2)))) {

		return false;
	}


	/* If the track is vertical */
	if (fabs(_phi - M_PI / 2) < 1E-10) {
		if (fabs(point->getX() - _start.getX()) < 1E-10)
			return true;
		else
			return false;
	}
	/* If the track is not vertical */
	else {
		m = sin(_phi) / cos(_phi);

		/* Use point-slope formula */
		if (fabs(point->getY() - (_start.getY() +
				m * (point->getX() - _start.getX()))) < 1e-10) {

			return true;
		}
		else
			return false;
	}
}


/**
 * Deletes each of this track's segments
 */
void Track::clearSegments() {
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
