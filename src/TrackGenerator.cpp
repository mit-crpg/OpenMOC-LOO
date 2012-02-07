/*
 * TrackGenerator.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "TrackGenerator.h"


/**
 * TrackGenerator constructor
 * @param geom a pointer to a geometry object
 * @param plotter a pointer to a plotting object
 * @param num_azim number of azimuthal angles
 * @param spacing track spacing
 */
TrackGenerator::TrackGenerator(Geometry* geom,
		const int num_azim, const double spacing, const int bitDim) {

	_geom = geom;
	_num_azim = num_azim/2.0;
	_spacing = spacing;

	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;

	_bit_length_x = int (bitDim*ratio) + 1;
	_bit_length_y = bitDim + 1;

	_pix_map_tracks = new int[_bit_length_x*_bit_length_y];
	_pix_map_segments = new int[_bit_length_x*_bit_length_y];

	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			_pix_map_segments[i * _bit_length_x + j] = -1;
		}
	}


	_x_pixel = double(_bit_length_x)/_width;
	_y_pixel = double(_bit_length_y)/_height;

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
}


/**
 * TrackGenerator destructor frees memory for all tracks
 */
TrackGenerator::~TrackGenerator() {
	delete [] _num_tracks;
	delete [] _azim_weights;
	delete _pix_map_segments;
	delete _pix_map_tracks;

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


/**
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

	/* Check to make sure that height, width of the geometry are nonzero */
	if (_geom->getHeight() <= 0 || _geom->getHeight() <= 0)
		log_printf(ERROR, "The total height and width of the geometry must be"
				"nonzero for track generation. Please specify the height and "
				"width in the geometry input file.");


	try {
		log_printf(NORMAL, "Computing azimuthal angles and track spacings...");

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


		/* Determine azimuthal angles and track spacing */
		for (int i = 0; i < _num_azim; i++) {

			/* desired angle */
			double phi = 2.0 * M_PI / iazim * (0.5 + i);

			/* num intersections with x,y-axes */
			_num_x[i] = (int) (fabs(width / _spacing * sin(phi))) + 1;
			_num_y[i] = (int) (fabs(height / _spacing * cos(phi))) + 1;

			/* total num of tracks */
			_num_tracks[i] = _num_x[i] + _num_y[i];

			/* effective/actual angle (not the angle we desire, but close */
			phi_eff[i] = atan((height * _num_x[i]) / (width * _num_y[i]));

			/* fix angles in range(pi/2, pi) */
			if (phi > M_PI / 2)
				phi_eff[i] = M_PI - phi_eff[i];

			/* Effective track spacing (not spacing we desire, but close */
			dx_eff[i] = (width / _num_x[i]);
			dy_eff[i] = (height / _num_y[i]);
			d_eff[i] = (dx_eff[i] * sin(phi_eff[i]));
		}

		/* Compute azimuthal weights */
		for (int i = 0; i < _num_azim; i++) {

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

		log_printf(NORMAL, "Generating track start and end points...");

		/* Compute track starting and end points */
		for (int i = 0; i < _num_azim; i++) {

			/* Tracks for azimuthal angle i */
			_tracks[i] = new Track[_num_tracks[i]];

			/* Compute start points for tracks starting on x-axis */
			for (int j = 0; j < _num_x[i]; j++)
				_tracks[i][j].getStart()->setCoords(dx_eff[i] * (0.5 + j), 0);

			/* Compute start points for tracks starting on y-axis */
			for (int j = 0; j < _num_y[i]; j++) {

				/* If track points to the upper right */
				if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) > 0)
					_tracks[i][_num_x[i]+j].getStart()->setCoords(0,
													dy_eff[i] * (0.5 + j));

				/* If track points to the upper left */
				else if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) < 0)
					_tracks[i][_num_x[i]+j].getStart()->setCoords(width,
														dy_eff[i] * (0.5 + j));
			}

			/* Compute the exit points for each track */
			for (int j = 0; j < _num_tracks[i]; j++) {

				/* Set the track's end point */
				Point* start = _tracks[i][j].getStart();
				Point* end = _tracks[i][j].getEnd();
				computeEndPoint(start, end, phi_eff[i], width, height);

				/* Set the track's azimuthal weight */
				_tracks[i][j].setWeight(_azim_weights[i]);
				_tracks[i][j].setPhi(phi_eff[i]);
			}
		}

		/* Recalibrate track start and end points to geometry global origin */
		//FIXME: This could be more efficiently done when start/end points are set
		for (int i = 0; i < _num_azim; i++) {
			for (int j = 0; j < _num_tracks[i]; j++) {
				double x0 = _tracks[i][j].getStart()->getX();
				double y0 = _tracks[i][j].getStart()->getY();
				double x1 = _tracks[i][j].getEnd()->getX();
				double y1 = _tracks[i][j].getEnd()->getY();
				double new_x0 = x0 - _geom->getWidth()/2.0;
				double new_y0 = y0 - _geom->getHeight()/2.0;
				double new_x1 = x1 - _geom->getWidth()/2.0;
				double new_y1 = y1 - _geom->getHeight()/2.0;
				double phi = _tracks[i][j].getPhi();
				_tracks[i][j].setValues(new_x0, new_y0, new_x1, new_y1, phi);

				/*
				 * Add line to _pix_map segments bitmap array.
				 * Note conversion from geometric coordinate system to
				 * bitmap coordinates.
				 */

				LineFct( new_x0*_x_pixel + _bit_length_x/2,
						-new_y0*_y_pixel + _bit_length_y/2,
						new_x1*_x_pixel + _bit_length_x/2,
						-new_y1*_y_pixel + _bit_length_y/2,
						_pix_map_tracks);
			}
		}

		return;
	}

	catch (std::exception &e) {
		log_printf(ERROR, "Unable to allocate memory needed to generate tracks"
				".Backtrace:\n%s", e.what());
	}

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
		const double phi, const double width, const double height) {

	double m = sin(phi) / cos(phi); 		/* slope */
	double yin = start->getY(); 			/* y-coord */
	double xin = start->getX(); 			/* x-coord */

	try {
		Point *points = new Point[4];

		/* Determine all possible points */
		points[0].setCoords(0, yin - m * xin);
		points[1].setCoords(width, yin + m * (width - xin));
		points[2].setCoords(xin - yin / m, 0);
		points[3].setCoords(xin - (yin - height) / m, height);

		/* For each of the possible intersection points */
		for (int i = 0; i < 4; i++) {
			/* neglect the trivial point (xin, yin) */
			if (points[i].getX() == xin && points[i].getY() == yin) { }

			/* The point to return will be within the bounds of the cell */
			else if (points[i].getX() >= 0 && points[i].getX() <= width
					&& points[i].getY() >= 0 && points[i].getY() <= height) {
				end->setCoords(points[i].getX(), points[i].getY());
			}
		}

		delete[] points;
		return;
	}

	catch (std::exception &e) {
		log_printf(ERROR, "Unable to allocate memory for intersection points "
				"in computeEndPoint method of TrackGenerator. Backtrace:\n%s",
				e.what());
	}

}


/**
 *  Implements reflective boundary conditions by setting the incoming
 *  and outgoing tracks for each track using a special indexing scheme
 *  into the trackgenerator's 2d jagged array of tracks
 */
void TrackGenerator::makeReflective() {
	log_printf(NORMAL, "Creating reflective boundary conditions...");

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
					curr[j].setReflIn(false);
					refl[j].setReflIn(false);

					curr[j].setTrackOut(&refl[2 * nxi - 1 - j]);
					refl[2 * nxi - 1 - j].setTrackIn(&curr[j]);
					curr[j].setReflOut(false);
					refl[2 * nxi - 1 - j].setReflIn(true);
				}

				/* Left hand side to right hand side */
				else if (j < nyi) {
					curr[j].setTrackIn(&refl[j - nxi]);
					refl[j - nxi].setTrackOut(&curr[j]);
					curr[j].setReflIn(true);
					refl[j - nxi].setReflOut(false);

					curr[j].setTrackOut(&refl[j + nxi]);
					refl[j + nxi].setTrackIn(&curr[j]);
					curr[j].setReflOut(false);
					refl[j + nxi].setReflIn(true);
				}

				/* Left hand side to top (j > ny) */
				else {
					curr[j].setTrackIn(&refl[j - nxi]);
					refl[j - nxi].setTrackOut(&curr[j]);
					curr[j].setReflIn(true);
					refl[j - nxi].setReflOut(false);

					curr[j].setTrackOut(&refl[2 * nti - nxi - j - 1]);
					refl[2 * nti - nxi - j - 1].setTrackOut(&curr[j]);
					curr[j].setReflOut(true);
					refl[2 * nti - nxi - j - 1].setReflOut(true);
				}
			}

			/* More tracks starting on y-axis than on x-axis */
			else {
				/* Bottom to top */
				if (j < nxi - nyi) {
					curr[j].setTrackIn(&refl[j]);
					refl[j].setTrackIn(&curr[j]);
					curr[j].setReflIn(false);
					refl[j].setReflIn(false);

					curr[j].setTrackOut(&refl[nti - (nxi - nyi) + j]);
					refl[nti - (nxi - nyi) + j].setTrackOut(&curr[j]);
					curr[j].setReflOut(true);
					refl[nti - (nxi - nyi) + j].setReflOut(true);
				}

				/* Bottom to right hand side */
				else if (j < nxi) {
					curr[j].setTrackIn(&refl[j]);
					refl[j].setTrackIn(&curr[j]);
					curr[j].setReflIn(false);
					refl[j].setReflIn(false);

					curr[j].setTrackOut(&refl[nxi + (nxi - j) - 1]);
					refl[nxi + (nxi - j) - 1].setTrackIn(&curr[j]);
					curr[j].setReflOut(false);
					refl[nxi + (nxi - j) - 1].setReflIn(true);
				}

				/* Left-hand side to top (j > nx) */
				else {
					curr[j].setTrackIn(&refl[j - nxi]);
					refl[j - nxi].setTrackOut(&curr[j]);
					curr[j].setReflIn(true);
					refl[j - nxi].setReflOut(false);

					curr[j].setTrackOut(&refl[nyi + (nti - j) - 1]);
					refl[nyi + (nti - j) - 1].setTrackOut(&curr[j]);
					curr[j].setReflOut(true);
					refl[nyi + (nti - j) - 1].setReflOut(true);
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

	/* Loop over all tracks */
	for (int i = 0; i < _num_azim; i++) {
		phi = _tracks[i][0].getPhi();
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		for (int j = 0; j < _num_tracks[i]; j++){
			_geom->segmentize(&_tracks[i][j]);
			plotSegmentsBitMap(&_tracks[i][j], sin_phi, cos_phi);
		}
	}

	return;
}


/**
 * Plot tracks in tiff file using fast drawing method
 */
void TrackGenerator::plotTracksTiff() {
	log_printf(NORMAL, "Writing tracks bitmap to tiff...");

	/* Create Magick image and open pixels for viewing/changing  */
	Magick::Image image_tracks(Magick::Geometry(_bit_length_x,_bit_length_y), "black");
	image_tracks.type(Magick::TrueColorType);
	Magick::Pixels view(image_tracks);

	/* Convert _pix_map_tracks bitmap array to Magick bitmap pixel array. */
	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			if (_pix_map_tracks[i * _bit_length_x + j] != 1){
				*(view.get(i,j,1,1)) = Magick::Color("white");
			}
		}
	}

	/* Close pixel viewing/changing */
	view.sync();

	/* Write pixel bitmap to tiff file */
	image_tracks.write("tracks.tiff");
}


/**
 * Plots track segments in a _pix_map_segments bitmap array on the fly
 */
void TrackGenerator::plotSegmentsBitMap(Track* track, double sin_phi, double cos_phi){

	/* Initialize variables */
	double start_x, start_y, end_x, end_y;
	int num_segments;

	/* Set first segment start point and get the number of tracks*/
	start_x = track->getStart()->getX();
	start_y = track->getStart()->getY();
	num_segments = track->getNumSegments();

	/* loop over segments and write to _pix_map_segments bitmap array */
	for (int k=0; k < num_segments; k++){
		end_x = start_x + cos_phi*track->getSegment(k)->_length;
		end_y = start_y + sin_phi*track->getSegment(k)->_length;

		/*
		 * Add line to _pix_map segments bitmap array.
		 * Note conversion from geometric coordinate system to
		 * bitmap coordinates.
		 */
		LineFct(start_x*_x_pixel + _bit_length_x/2,
				-start_y*_y_pixel + _bit_length_y/2,
				end_x*_x_pixel + _bit_length_x/2,
				-end_y*_y_pixel + _bit_length_y/2,
				_pix_map_segments,
				track->getSegment(k)->_region_id);

		start_x = end_x;
		start_y = end_y;
	}
}



/**
 * Plot segments in tiff file
 */
void TrackGenerator::plotSegmentsTiff(){
	log_printf(NORMAL, "Writing segments bitmap to tiff...");

	/* Create Magick image and open pixels for viewing/changing */
	Magick::Image image_segments(Magick::Geometry(_bit_length_x,_bit_length_y), "white");
	image_segments.type(Magick::TrueColorType);
	Magick::Pixels view(image_segments);

	/*
	 * Convert _pix_map_segments bitmap array to Magick bitmap pixel
	 * color array
	 */
	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			switch (_pix_map_segments[i * _bit_length_x + j]){
			case 0:
				*(view.get(i,j,1,1)) = Magick::Color("indigo");
				break;
			case 1:
				*(view.get(i,j,1,1)) = Magick::Color("red");
				break;
			case 2:
				*(view.get(i,j,1,1)) = Magick::Color("blue");
				break;
			case 3:
				*(view.get(i,j,1,1)) = Magick::Color("green");
				break;
			case 4:
				*(view.get(i,j,1,1)) = Magick::Color("magenta");
				break;
			case 5:
				*(view.get(i,j,1,1)) = Magick::Color("orange");
				break;
			case 6:
				*(view.get(i,j,1,1)) = Magick::Color("maroon");
				break;
			case 7:
				*(view.get(i,j,1,1)) = Magick::Color("orchid");
				break;
			case 8:
				*(view.get(i,j,1,1)) = Magick::Color("blue violet");
				break;
			case 9:
				*(view.get(i,j,1,1)) = Magick::Color("crimson");
				break;
			case 10:
				*(view.get(i,j,1,1)) = Magick::Color("salmon");
				break;
			case 11:
				*(view.get(i,j,1,1)) = Magick::Color("gold");
				break;
			case 12:
				*(view.get(i,j,1,1)) = Magick::Color("blue violet");
				break;
			case 13:
				*(view.get(i,j,1,1)) = Magick::Color("orange red");
				break;
			case 14:
				*(view.get(i,j,1,1)) = Magick::Color("spring green");
				break;
			case 15:
				*(view.get(i,j,1,1)) = Magick::Color("pale green");
				break;
			}
		}
	}

	/* close pixel viewing/changing */
	view.sync();

	/* write Magick pixel color array to tiff file */
	image_segments.write("segments.tiff");
}


/**
 *  finds the sign of a value
 *  @param a find the sign of a
 */
int TrackGenerator::sgn (long a) {
	if (a > 0) return +1;
	else if (a < 0) return -1;
	else return 0;
}

/**
 * Line drawing algorithm. Takes in the start (a,b) and end (c,d) coordinates
 * of line, pointer to _pix_map bitmap array (pixMap), and line color (col).
 * "Draws" the line on _pix_map bitmap array.
 * Algorithm taken from http://www.cprogramming.com/tutorial/tut3.html
*/
void TrackGenerator::LineFct(int a, int b, int c, int d, int* pixMap, int col) {
	long u,s,v,d1x,d1y,d2x,d2y,m,n;
	int  i;
	u   = c-a;
	v   = d-b;
	d1x = sgn(u);
	d1y = sgn(v);
	d2x = sgn(u);
	d2y = 0;
	m   = abs(u);
	n   = abs(v);
	if (m<=n) {
		d2x = 0;
		d2y = sgn(v);
		m   = abs(v);
	    n   = abs(u);
	}
	s = (int)(m / 2);
	for (i=0;i<round(m);i++) {
		pixMap[a * _bit_length_x + b] = col;
		s += n;
		if (s >= m) {
			s -= m;
			a += d1x;
			b += d1y;
		}
		else {
			a += d2x;
			b += d2y;
		}
	}
}


