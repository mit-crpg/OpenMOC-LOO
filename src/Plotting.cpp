/*
 * Plotting.cpp
 *
 *  Created on: Feb 1, 2012
 *      Author: samuelshaner
 */

#include "Plotting.h"


Plotting::Plotting(Geometry* geom){

	_geom = geom;

	/* get width and height */
	double width = _geom->getWidth();
	double height = _geom->getHeight();
	double ratio = width/height;

	_bit_length_x = int (1000*ratio);
	_bit_length_y = 1000;

	_x_pixel = _bit_length_x/width;
	_y_pixel = _bit_length_y/height;
}

Plotting::~Plotting(){
}


void Plotting::plotSegments(TrackGenerator* track_generator){
	log_printf(NORMAL, "Creating tiff plot of segments...");

	std::list<Magick::Drawable> drawListRed0;
	std::list<Magick::Drawable> drawListBlue1;
	std::list<Magick::Drawable> drawListGreen2;
	std::list<Magick::Drawable> drawListPink3;
	std::list<Magick::Drawable> drawListOrange4;
	std::list<Magick::Drawable> drawListMaroon5;
	std::list<Magick::Drawable> drawListOrchid6;
	std::list<Magick::Drawable>* _draw_lists[7] = {&drawListRed0, &drawListBlue1,
			&drawListGreen2, &drawListPink3, &drawListOrange4, &drawListMaroon5, &drawListOrchid6};
	_draw_lists[0]->push_back(Magick::DrawableStrokeColor("Red"));
	_draw_lists[0]->push_back(Magick::DrawableFillColor("Red"));
	_draw_lists[1]->push_back(Magick::DrawableStrokeColor("Blue"));
	_draw_lists[1]->push_back(Magick::DrawableFillColor("Blue"));
	_draw_lists[2]->push_back(Magick::DrawableStrokeColor("Green"));
	_draw_lists[2]->push_back(Magick::DrawableFillColor("Green"));
	_draw_lists[3]->push_back(Magick::DrawableStrokeColor("Pink"));
	_draw_lists[3]->push_back(Magick::DrawableFillColor("Pink"));
	_draw_lists[4]->push_back(Magick::DrawableStrokeColor("Orange"));
	_draw_lists[4]->push_back(Magick::DrawableFillColor("Orange"));
	_draw_lists[5]->push_back(Magick::DrawableStrokeColor("Maroon"));
	_draw_lists[5]->push_back(Magick::DrawableFillColor("Maroon"));
	_draw_lists[6]->push_back(Magick::DrawableStrokeColor("Orchid"));
	_draw_lists[6]->push_back(Magick::DrawableFillColor("Orchid"));

	for (int i = 0; i < 7; i++){
		_draw_lists[i]->push_back(Magick::DrawableTranslation(_bit_length_x/2,_bit_length_y/2));
		_draw_lists[i]->push_back(Magick::DrawableStrokeWidth(0));
		_draw_lists[i]->push_back(Magick::DrawableStrokeAntialias(false));
	}


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
			log_printf(DEBUG, "Now printing segments for track i: %d, j: %d", i, j);
			for (int k=0; k < num_segments; k++){
				 end_x = start_x + cos_phi*tracks[i][j].getSegment(k)->_length;
				 end_y = start_y + sin_phi*tracks[i][j].getSegment(k)->_length;
				 log_printf(DEBUG, "start_x: %f, start_y: %f, end_x: %f, end_y: %f, region: %d",
						 start_x, start_y, end_x, end_y, tracks[i][j].getSegment(k)->_region_id);
//				 log_printf(DEBUG, "start_x: %f, start_y: %f, end_x: %f, end_y: %f, region: %d",
//				 						 start_x*_x_pixel, start_y*_y_pixel, end_x*_x_pixel, end_y*_y_pixel,
//				 						 tracks[i][j].getSegment(k)->_region_id);
				 _draw_lists[tracks[i][j].getSegment(k)->_region_id % 7]->push_back(Magick::DrawableLine
						 (start_x*_x_pixel, -start_y*_y_pixel, end_x*_x_pixel, -end_y*_y_pixel));

				 start_x = end_x;
				 start_y = end_y;
			}
		}
	}

	/* make image and write tiff file */
	Magick::Image image_segments(Magick::Geometry(_bit_length_x,_bit_length_y), Magick::Color("white"));

	for (int i = 0; i < 7; i++){
		image_segments.draw(*_draw_lists[i]);
	}

	image_segments.write("segments.tiff");

}


/*
 * plot tracks in tiff file
 */
void Plotting::plotTracksTiff(TrackGenerator* track_generator) {
	log_printf(NORMAL, "Creating tiff plot of tracks...");

	// initialize image
	Magick::Image image_tracks(Magick::Geometry(_bit_length_x,_bit_length_y), Magick::Color("white"));
	std::list<Magick::Drawable> drawList;

	// translate (width/2,height/2), rotate -90
	drawList.push_back(Magick::DrawableTranslation(_bit_length_x/2,_bit_length_y/2));

	drawList.push_back(Magick::DrawableStrokeColor("black"));
	drawList.push_back(Magick::DrawableFillColor("black"));
	drawList.push_back(Magick::DrawableStrokeWidth(0));
	drawList.push_back(Magick::DrawableStrokeAntialias(false));

	Track** tracks = track_generator->getTracks();
	int* num_tracks = track_generator->getNumTracks();

	// append tracks to drawList
	for (int i=0; i < track_generator->getNumAzim(); i++) {
		for (int j=0; j < num_tracks[i]; j++) {
			drawList.push_back(Magick::DrawableLine( tracks[i][j].getStart()->getX()*_x_pixel,
					-tracks[i][j].getStart()->getY()*_y_pixel
					, tracks[i][j].getEnd()->getX()*_x_pixel,
					-tracks[i][j].getEnd()->getY()*_y_pixel));
		}
	}

	// output tracks to tracks.tiff
	image_tracks.draw(drawList);
	image_tracks.write("tracks.tiff");
}




