/*
 * Plotting.cpp
 *
 *  Created on: Feb 1, 2012
 *      Author: samuelshaner
 */

#include "Plotting.h"


Plotting::Plotting(Geometry* geom, int bitDim){

	_geom = geom;

	/* get width and height */
	double width = _geom->getWidth();
	double height = _geom->getHeight();
	double ratio = width/height;

	_bit_length_x = int (bitDim*ratio) + 1;
	_bit_length_y = bitDim + 1;

	_pix_map_tracks = new int[_bit_length_x*_bit_length_y];
	_pix_map_segments = new int[_bit_length_x*_bit_length_y];


	_x_pixel = double(_bit_length_x)/width;
	_y_pixel = double(_bit_length_y)/height;
}

Plotting::~Plotting(){
	delete _pix_map_tracks;
	delete _pix_map_segments;
}


/* plot segments in tiff file using fast drawing method */
void Plotting::plotSegments(TrackGenerator* track_generator){
	log_printf(NORMAL, "Creating tiff plot of segments...");

	/* generate _pix_map */

	/* loop over azim, tracks, segments */

	Track** tracks = track_generator->getTracks();
	int* num_tracks = track_generator->getNumTracks();
	double phi;

	double sin_phi;
	double cos_phi;
	double start_x, start_y, end_x, end_y;
	int num_segments;

	// append tracks to drawList
	for (int i=0; i < track_generator->getNumAzim(); i++) {
		for (int j=0; j < num_tracks[i]; j++) {
			phi = tracks[i][j].getPhi();
			sin_phi = sin(phi);
			cos_phi = cos(phi);
			start_x = tracks[i][j].getStart()->getX();
			start_y = tracks[i][j].getStart()->getY();
			num_segments = tracks[i][j].getNumSegments();

			for (int k=0; k < num_segments; k++){
				end_x = start_x + cos_phi*tracks[i][j].getSegment(k)->_length;
				end_y = start_y + sin_phi*tracks[i][j].getSegment(k)->_length;
				SegFct(start_x*_x_pixel + _bit_length_x/2,
						-start_y*_y_pixel + _bit_length_y/2,
						end_x*_x_pixel + _bit_length_x/2,
						-end_y*_y_pixel + _bit_length_y/2,
						tracks[i][j].getSegment(k)->_region_id);

				start_x = end_x;
				start_y = end_y;
			}
		}
	}


	for (int i=0; i < track_generator->getNumAzim(); i++) {
		for (int j=0; j < num_tracks[i]; j++) {
			LineFct( tracks[i][j].getStart()->getX()*_x_pixel + _bit_length_x/2,
					-tracks[i][j].getStart()->getY()*_y_pixel + _bit_length_y/2
					, tracks[i][j].getEnd()->getX()*_x_pixel + _bit_length_x/2,
					-tracks[i][j].getEnd()->getY()*_y_pixel + _bit_length_y/2, 1);
		}
	}

	/* write _pix_map to tiff file  */
	Magick::Image image_segments(Magick::Geometry(_bit_length_x,_bit_length_y), "white");

	image_segments.type(Magick::TrueColorType);

	Magick::Pixels view(image_segments);

	log_printf(NORMAL, "coloring pixels...");

	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			switch (_pix_map_segments[i * _bit_length_x + j]){
			case 0:
				break;
			case 1:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("red");
				break;
			case 2:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("blue");
				break;
			case 3:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("green");
				break;
			case 4:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("magenta");
				break;
			case 5:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("orange");
				break;
			case 6:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("maroon");
				break;
			case 7:
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("orchid");
				break;
			}
		}
	}

	view.sync();

	image_segments.write("segments.tiff");
}


/*
 * plot tracks in tiff file using fast drawing method
 */
void Plotting::plotTracksTiff(TrackGenerator* track_generator) {
	log_printf(NORMAL, "Creating tiff plot of tracks...");

	/* generate _pix_map */

	Track** tracks = track_generator->getTracks();
	int* num_tracks = track_generator->getNumTracks();

	// append tracks to drawList
	for (int i=0; i < track_generator->getNumAzim(); i++) {
		for (int j=0; j < num_tracks[i]; j++) {
			LineFct( tracks[i][j].getStart()->getX()*_x_pixel + _bit_length_x/2,
					-tracks[i][j].getStart()->getY()*_y_pixel + _bit_length_y/2
					, tracks[i][j].getEnd()->getX()*_x_pixel + _bit_length_x/2,
					-tracks[i][j].getEnd()->getY()*_y_pixel + _bit_length_y/2, 1);
		}
	}

	/* write _pix_map to tiff file  */
	Magick::Image image_tracks(Magick::Geometry(_bit_length_x,_bit_length_y), "white");

	image_tracks.type(Magick::TrueColorType);

	Magick::Pixels view(image_tracks);

	log_printf(NORMAL, "coloring pixels...");

	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			if (_pix_map_tracks[i * _bit_length_x + j] != 0){
				*(view.get(i,_bit_length_y - 1 - j,1,1)) = Magick::Color("black");
			}
		}
	}

	view.sync();

	image_tracks.write("tracks.tiff");
}

/* finds the sign of a value */
int Plotting::sgn (long a) {
	if (a > 0) return +1;
	else if (a < 0) return -1;
	else return 0;
}

/*
 * track drawing algorithm
 * taken from http://www.cprogramming.com/tutorial/tut3.html
*/
void Plotting::LineFct(int a, int b, int c, int d, int col) {
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
		_pix_map_tracks[a * _bit_length_x + b] = col;
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

/*
 * segment drawing algorithm
 * taken from http://www.cprogramming.com/tutorial/tut3.html
 */
void Plotting::SegFct(int a, int b, int c, int d, int col) {
	col = col % 7 + 1;
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
		_pix_map_segments[a * _bit_length_x + b] = col;
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



