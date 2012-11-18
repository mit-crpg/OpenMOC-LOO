/*
 * Mesh.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "Mesh.h"

Mesh::Mesh(){
}

Mesh::~Mesh(){
	delete [] _cells;
}

void Mesh::makeMeshCells(){

	/* make mesh cells */
	_cells = new MeshCell[_cell_width * _cell_height];

}

int Mesh::getCellWidth(){
	return _cell_width;
}

int Mesh::getCellHeight(){
	return _cell_height;
}

void Mesh::setCellWidth(int cellWidth){
	_cell_width = cellWidth;
}

void Mesh::setCellHeight(int cellHeight){
	_cell_height = cellHeight;
}

MeshCell* Mesh::getCells(int cell){
	return &_cells[cell];
}


void Mesh::setMultigroup(bool multigroup){
	_multigroup = multigroup;
}

bool Mesh::getMultigroup(){
	return _multigroup;
}

void Mesh::setPrintMatrices(bool printMatrices){
	_print_matrices = printMatrices;
}

bool Mesh::getPrintMatrices(){
	return _print_matrices;
}

/* given an x,y coordinate, find what mesh cell the point is in */
int Mesh::findMeshCell(double pointX, double pointY){

	double left;
	double top = _height / 2.0;
	bool flag = false;
	int cell = 0;

	for (int y = 0; y < _cell_height; y++){
		left = - _width / 2.0;
		cell = y * _cell_width;
		if (pointY <= top && pointY >= top - getCells(cell)->getHeight()){
			for (int x = 0; x < _cell_width; x++){
				cell = y * _cell_width + x;
				if (pointX >= left && pointX <= left + getCells(cell)->getWidth()){
					flag = true;
					break;
				}
				left = left + getCells(cell)->getWidth();
			}
			if (flag == true){
				break;
			}
		}
		top = top - getCells(cell)->getHeight();
	}

	return cell;
}

void Mesh::setWidth(double width){
	_width = width;
}

void Mesh::setHeight(double height){
	_height = height;
}

double Mesh::getWidth(){
	return _width;
}

double Mesh::getHeight(){
	return _height;
}

void Mesh::setCellBounds(){

	double x = -_width / 2.0;
	double y = _height / 2.0;

	/* loop over MeshCells and set bounds */
	for (int i = 0; i < _cell_height; i++){
		x = -_width / 2.0;
		y = y - _cells[i * _cell_width].getHeight();
		for (int j = 0; j < _cell_width; j++){
			_cells[i * _cell_width + j].setBounds(x,y);
			x = x + _cells[i * _cell_width + j].getWidth();
		}
	}

}

void Mesh::setFSRBounds(){

	int min;
	int max;
	int fsr;

	for (int i = 0; i < _cell_height * _cell_width; i++){
		min = _cells[i].getFSRs()->front();
		max = _cells[i].getFSRs()->front();

		std::vector<int>::iterator iter;
		for (iter = _cells[i].getFSRs()->begin(); iter != _cells[i].getFSRs()->end(); ++iter) {
			fsr = *iter;
			min = std::min(fsr, min);
			max = std::max(fsr, max);
		}

		_cells[i].setFSRStart(min);
		_cells[i].setFSREnd(max);
	}
}


/* Using an fsr_id and coordinate, find which surface a coordinate is on */
MeshSurface* Mesh::findMeshSurface(int fsr_id, LocalCoords* coord){
	MeshSurface* meshSurface = NULL;

	/* find which MeshCell fsr_id is in -> get meshSuface that coord is on*/
	for (int i = 0; i < _cell_width * _cell_height; i++){
		if (fsr_id >= _cells[i].getFSRStart() && fsr_id <= _cells[i].getFSREnd()){
			meshSurface = _cells[i].findSurface(coord, i);
			break;
		}

	}

	return meshSurface;
}

void Mesh::printBounds(){

	double* bounds;

	for (int i = 0; i < _cell_width * _cell_height; i++){
		bounds = _cells[i].getBounds();
		log_printf(NORMAL, "cell: %i bounds [%f, %f, %f, %f]", i, bounds[0], bounds[1], bounds[2], bounds[3]);
	}

}


void Mesh::printCurrents(){

	double current;

	for (int i = 0; i < _cell_width * _cell_height; i++){
		for (int surface = 0; surface < 8; surface++){
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
				current = _cells[i].getMeshSurfaces(surface)->getCurrent(group);
				log_printf(NORMAL, "cell: %i, surface: %i, group: %i, current: %f", i, surface, group, current);
			}
		}
	}
}


void Mesh::splitCorners(){

	log_printf(NORMAL, "Splitting corners...");

	MeshSurface* surfaceSide;
	MeshSurface* surfaceCorner1;

	int cell_width = getCellWidth();
	int cell_height = getCellHeight();
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	MeshSurface* surfaceSideNext;

	for (int x = 0; x < cell_width; x++){
		for (int y = 0; y < cell_height; y++){

			/* corner 4 */
			if (x > 0 && y < cell_height - 1){
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(1);
				surfaceCorner1 = meshCell->getMeshSurfaces(4);
				meshCellNext = &_cells[(y+1)*cell_width + x];
				surfaceSideNext = meshCellNext->getMeshSurfaces(0);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}
			else{
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(1);
				surfaceCorner1 = meshCell->getMeshSurfaces(4);
				surfaceSideNext = meshCell->getMeshSurfaces(0);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}


			/* corner 5 */
			if (x < cell_width - 1 && y < cell_height - 1){
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(1);
				surfaceCorner1 = meshCell->getMeshSurfaces(5);
				meshCellNext = &_cells[(y+1)*cell_width + x];
				surfaceSideNext = meshCellNext->getMeshSurfaces(2);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}
			else{
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(1);
				surfaceCorner1 = meshCell->getMeshSurfaces(5);
				surfaceSideNext = meshCell->getMeshSurfaces(2);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}


			/* corner 6 */
			if (x < cell_width - 1 && y > 0){
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(2);
				surfaceCorner1 = meshCell->getMeshSurfaces(6);
				meshCellNext = &_cells[y*cell_width + x + 1];
				surfaceSideNext = meshCellNext->getMeshSurfaces(3);

					for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
						surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
						surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					}
			}
			else{
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(2);
				surfaceCorner1 = meshCell->getMeshSurfaces(6);
				surfaceSideNext = meshCell->getMeshSurfaces(3);

					for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
						surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
						surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					}
			}


			/* corner 7 */
			if (x > 0 && y > 0){
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(0);
				surfaceCorner1 = meshCell->getMeshSurfaces(7);
				meshCellNext = &_cells[y*cell_width + x - 1];
				surfaceSideNext = meshCellNext->getMeshSurfaces(3);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}
			else{
				meshCell = &_cells[y*cell_width + x];
				surfaceSide = meshCell->getMeshSurfaces(0);
				surfaceCorner1 = meshCell->getMeshSurfaces(7);
				surfaceSideNext = meshCell->getMeshSurfaces(3);

				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					surfaceSide->incrementCurrent(surfaceCorner1->getCurrent(group), group);
					surfaceSideNext->incrementCurrent(surfaceCorner1->getCurrent(group), group);
				}
			}
		}
	}
}

void Mesh::computeTotCurrents(){

	MeshCell* meshCell;
	MeshSurface* surfaceSide;
	double sum_cur;

	int cell_width = getCellWidth();
	int cell_height = getCellHeight();

	/* loop over cells */
	for (int i = 0; i < cell_width*cell_height; i++){

		/* get mesh cell */
		meshCell = &_cells[i];

		/* loop over surfaces*/
		for (int i = 0; i < 4; i++){

			/* get mesh surface */
			surfaceSide = meshCell->getMeshSurfaces(i);

			/* set current tally to 0 */
			sum_cur = 0.0;

			/* loop over energy groups */
			for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
				sum_cur += surfaceSide->getCurrent(e);
			}

			/* set current at group 0 to total current */
			surfaceSide->setCurrent(sum_cur,0);

		}
	}
}


