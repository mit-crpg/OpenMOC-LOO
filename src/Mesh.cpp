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

    /* set mesh surface cell id's */
    for (int i = 0; i < _cell_width * _cell_height; i++){
        for (int s = 0; s < 8; s++){
            _cells[i].getMeshSurfaces(s)->setCellId(i);
        }
    }

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


int Mesh::findMeshCell(double pointX, double pointY)
{
    double left = - _width / 2.0;
    double bottom = - _height / 2.0;
    double height, width;

    for (int y = 0; y < _cell_height; y++){
        height = getCells(y * _cell_width)->getHeight();

        if ((pointY > bottom - ON_LATTICE_CELL_THRESH) && 
            (pointY < bottom + height + ON_LATTICE_CELL_THRESH))
        {
            for (int x = 0; x < _cell_width; x++)
            {
                width = getCells(y * _cell_width + x)->getWidth();

                if ((pointX > left - ON_LATTICE_CELL_THRESH)
                    && (pointX < left + width + ON_LATTICE_CELL_THRESH))
                    return y * _cell_width + x;
                left += width;
            }
        }
        bottom += height;
    }

    return -1;
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

void Mesh::setFSRBounds(boundaryType left, boundaryType right, boundaryType bottom, boundaryType top){

    _fsr_indices = new int        [2 * _cell_width * _cell_height];
    _cell_bounds = new double     [4* _cell_width * _cell_height];
    _surfaces =  new MeshSurface *[8 * _cell_width * _cell_height];

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
        _fsr_indices[2*i] = min;
        _fsr_indices[2*i+1] = max;

        for (int iii = 0; iii < 4; iii++)
            _cell_bounds[i * 4 + iii] = _cells[i].getBounds()[iii];

        for (int iii = 0; iii < 8; iii++)
            _surfaces[i * 8 + iii] = _cells[i].getMeshSurfaces(iii);
    }

    for (int x = 0; x < _cell_width; x++){
        for (int y = 0; y < _cell_height; y++){

            /* left */
            if (x == 0){
                if (left == VACUUM){
                    _cells[y*_cell_width+x].getMeshSurfaces(0)->setBoundary(VACUUM);
                }
                else{
                    _cells[y*_cell_width+x].getMeshSurfaces(0)->setBoundary(REFLECTIVE);
                }
            }

            /* right */
            if (x == _cell_width-1){
                if (right == VACUUM){
                    _cells[y*_cell_width+x].getMeshSurfaces(2)->setBoundary(VACUUM);
                }
                else{
                    _cells[y*_cell_width+x].getMeshSurfaces(2)->setBoundary(REFLECTIVE);
                }
            }

            /* bottom */
            if (y == _cell_height - 1){
                if (bottom == VACUUM){
                    _cells[y*_cell_width+x].getMeshSurfaces(1)->setBoundary(VACUUM);
                }
                else{
                    _cells[y*_cell_width+x].getMeshSurfaces(1)->setBoundary(REFLECTIVE);
                }
            }

            /* top */
            if (y == 0){
                if (top == VACUUM){
                    _cells[y*_cell_width+x].getMeshSurfaces(3)->setBoundary(VACUUM);
                }
                else{
                    _cells[y*_cell_width+x].getMeshSurfaces(3)->setBoundary(REFLECTIVE);
                }
            }

        }
    }


}


/* Using an fsr_id and coordinate, find which surface a coordinate is on */
int Mesh::findMeshSurface(int fsr_id, LocalCoords* coord)
{
    int surface = -1;
    double x = coord->getX();
    double y = coord->getY();

    double limit = 1e-10;

    /* find which MeshCell fsr_id is in -> get meshSuface that coord is on*/
    for (int i = 0; i < _cell_width * _cell_height; i++)
    {
        if (fsr_id >= _fsr_indices[2*i] && fsr_id <= _fsr_indices[2*i+1])
        {
            /* find which surface coord is on */
            /* left */
            if (fabs(x - _cell_bounds[i*4+0]) < limit)
            {
                if (fabs(y - _cell_bounds[i*4+1]) > limit && 
                    fabs(y - _cell_bounds[i*4+3]) > limit)
                    surface = i*8+0;
                else if (fabs(y - _cell_bounds[i*4+3]) < limit)
                    surface = i*8+7;
                else
                    surface = i*8+4;
            }
            /* right */
            else if (fabs(x - _cell_bounds[i*4+2]) < limit)
            {
                if (fabs(y - _cell_bounds[i*4+1]) > limit && 
                    fabs(y - _cell_bounds[i*4+3]) > limit)
                    surface = i*8+2;
                else if (fabs(y - _cell_bounds[i*4+3]) < limit)
                    surface = i*8+6;
                else
                    surface = i*8+5;
            }
            /* top */
            else if (fabs(y - _cell_bounds[i*4+3]) < limit)
                surface = i*8+3;
            /* bottom */
            else if (fabs(y - _cell_bounds[i*4+1]) < limit)
                surface = i*8+1;

            break;
        }
    }

    return surface;
}

void Mesh::printBounds()
{
    double* bounds;

    for (int i = 0; i < _cell_width * _cell_height; i++)
    {
        bounds = _cells[i].getBounds();
        log_printf(NORMAL, "cell: %i bounds [%f, %f, %f, %f]", 
                   i, bounds[0], bounds[1], bounds[2], bounds[3]);
    }
    
    return;
}


void Mesh::printCurrents()
{
    double current;

    for (int i = 0; i < _cell_width * _cell_height; i++)
    {
        for (int surface = 0; surface < 8; surface++)
        {
            for (int group = 0; group < NUM_ENERGY_GROUPS; group++)
            {
                current = _cells[i].getMeshSurfaces(surface)->getCurrent(group);
                log_printf(NORMAL, 
                           "cell: %i, surface: %i, group: %i, current: %f", 
                           i, surface, group, current);
            }
        }
    }
}


void Mesh::setBoundary(boundaryType boundary, int s){
    assert(s >= 0);
    assert(s < 4);
    _boundary[s] = boundary;
}


boundaryType Mesh::getBoundary(int s){
    assert(s >= 0);
    assert(s < 4);
    return _boundary[s];
}

void Mesh::splitCornerCurrents()
{
    MeshSurface* surfaceCorner;
    MeshSurface* surfaceSide;
    MeshSurface* surfaceSideNext;

    int cw = getCellWidth();
    int ch = getCellHeight();
    MeshCell* meshCell;
    MeshCell* meshCellNext;
    double currents[NUM_ENERGY_GROUPS];
    double f = 1.0; 

#if 0
    int min_x[] = {0, -1, -1, 0};
    int max_x[] = {cw, cw-1, cw-1, cw};
    int min_y[] = {-1, -1, 0, 0};
    int max_y[] = {ch-1, ch-1, ch, ch};
#else
    int min_x[] = {-1,   -1,   -1,   0};
    int max_x[] = {cw,   cw,   cw-1, cw};
    int min_y[] = {-1,   -1,   -1,    -1};
    int max_y[] = {ch-1, ch-1, ch,   ch};
#endif

    int next_x[] = {0, 0, 1, -1};
    int next_y[] = {1, 1, 0, 0};

    int surf[]      = {1, 1, 2, 0};
    int crn[]       = {4, 5, 6, 7};
    int next_surf[] = {0, 2, 3, 3};

    for (int x = 0; x < cw; x++){
        for (int y = 0; y < ch; y++){
            meshCell = &_cells[y * cw + x];

            log_printf(DEBUG, "cell %d's four corners' have current: %e %e"
                       " %e %e", y * cw + x, 
                       meshCell->getMeshSurfaces(4)->getCurrent(0),
                       meshCell->getMeshSurfaces(5)->getCurrent(0), 
                       meshCell->getMeshSurfaces(6)->getCurrent(0),
                       meshCell->getMeshSurfaces(7)->getCurrent(0));
					   
            for (int i = 0; i < 4; i++)
            {
                if (x > min_x[i] && x < max_x[i] && 
                    y > min_y[i] && y < max_y[i])
                {
                    meshCellNext = 
                        &_cells[(y + next_y[i]) * cw + x + next_x[i]];
                }
                else 
                    meshCellNext = meshCell;

                surfaceSide = meshCell->getMeshSurfaces(surf[i]);
                surfaceSideNext = meshCellNext->getMeshSurfaces(next_surf[i]);
                surfaceCorner = meshCell->getMeshSurfaces(crn[i]);
                for (int group = 0; group < NUM_ENERGY_GROUPS; group++)
                    currents[group] = f * surfaceCorner->getCurrent(group);

                surfaceSide->incrementCurrents(currents);
                surfaceSideNext->incrementCurrents(currents); 			
            }
        }
    }
    return;
}
*/

void Mesh::splitCornerQuadWeights()
{
    MeshSurface* surfaceCorner;
    MeshSurface *surface1, *surface2, *surface_next1, *surface_next2;

    int cw = getCellWidth();
    int ch = getCellHeight();
    MeshCell *meshCell, *meshCellNext, *meshCellNext2; 
    double current;
    double f = 0.50; 

    /* associated with this current cell */
    int surf[]      = {0, 2, 2, 0};
    int surf2[]     = {1, 1, 3, 3};

    /* associated with cells in the x-plane */
    int min_x[] = {0,    -1,   -1,   0};
    int max_x[] = {cw,   cw-1, cw-1, cw};
    int min_y[] = {-1,   -1,   0,    0};
    int max_y[] = {ch-1, ch-1, ch,   ch};

    int next_x[]    = {-1, 1, 1, -1};
    int next_surf[] = {1,  1, 3, 3};
    
    int next_y2[]   = {1, 1, -1, -1};
    int next_surf2[]= {0, 2, 2, 0};

    /* number of quadrature currents that we care to split.*/
    int nq = 2; 

    /* j is the index of the quadrature, counter_j stores the
     * reflective quadrature index */
    int counter_j[] = {1, 0, 3, 2};

    for (int y = 0; y < ch; y++)
    {
        for (int x = 0; x < cw; x++)
        {
            int ii = y * cw + x;
            meshCell = &_cells[ii];

            for (int i = 0; i < 4; i++)
            {
                surfaceCorner = meshCell->getMeshSurfaces(i + 4);

                /* perform splitting inside of this cell: distribute
                 * evenly to the two surfaces surrounding the
                 * corner */
                surface1 = meshCell->getMeshSurfaces(surf[i]);
                surface2 = meshCell->getMeshSurfaces(surf2[i]);

                for (int j = 0; j < nq; j++)
                {
                    current = f * surfaceCorner->getTotalWt(j);
                    surface1->incrementTotalWt(current, j);
                    surface2->incrementTotalWt(current, j);
                }

                /* effects on the neighboring cells */
                /* if left or right cell exists */
                if ((x > min_x[i]) && (x < max_x[i]))
                {
                    meshCellNext = 
                        &_cells[ii + next_x[i]];
                    surface_next1 = meshCellNext->getMeshSurfaces(next_surf[i]);
                   
                     for (int j = 0; j < nq; j++)
                     {
                         current = f * surfaceCorner->getTotalWt(j);
                         surface_next1->incrementTotalWt(current, j);
                     }
                } 
                /* only need to worry about this when reflective */
                else if (((x == min_x[i]) && (_boundary[0] == REFLECTIVE)) 
                         || ((x == max_x[i]) && (_boundary[2] == REFLECTIVE)))
                {
                    surface_next1 = meshCell->getMeshSurfaces(next_surf[i]);
                   
                     for (int j = 0; j < nq; j++)
                     {
                         current = f * surfaceCorner->getTotalWt(j);
                             surface_next1->incrementTotalWt
                                 (current, counter_j[j]);
                     }
                }



                if ((y > min_y[i]) && (y < max_y[i])) 
                {
                     meshCellNext2 = 
                        &_cells[ii + next_y2[i] * cw];
                     surface_next2 = meshCellNext2->getMeshSurfaces
                        (next_surf2[i]);

                     for (int j = 0; j < nq; j++)
                     {
                         current = f * surfaceCorner->getTotalWt(j);
                         surface_next2->incrementTotalWt(current, j);
                     }
                } 
                else if (((y == min_y[i]) && (_boundary[3] == REFLECTIVE)) 
                         || ((y == max_y[i]) && (_boundary[1] == REFLECTIVE)))
                {
                     surface_next2 = meshCell->getMeshSurfaces(next_surf2[i]);

                     for (int j = 0; j < nq; j++)
                     {
                         current = f * surfaceCorner->getTotalWt(j);
                         surface_next2->incrementTotalWt
                             (current, counter_j[j]);
                     }
                } 



            }
        }
    }
    return;
}

void Mesh::splitCornerQuadCurrents()
{
    MeshSurface* surfaceCorner;
    MeshSurface *surface1, *surface2, *surface_next1, *surface_next2;

    int cw = getCellWidth();
    int ch = getCellHeight();
    MeshCell *meshCell, *meshCellNext, *meshCellNext2; 
    double current;
    double f = 0.50; 

    /* associated with this current cell */
    int surf[]      = {0, 2, 2, 0};
    int surf2[]     = {1, 1, 3, 3};

    /* associated with cells in the x-plane */
    int min_x[] = {0,    -1,   -1,   0};
    int max_x[] = {cw,   cw-1, cw-1, cw};
    int min_y[] = {-1,   -1,   0,    0};
    int max_y[] = {ch-1, ch-1, ch,   ch};

    int next_x[]    = {-1, 1, 1, -1};
    int next_surf[] = {1,  1, 3, 3};
    
    int next_y2[]   = {1, 1, -1, -1};
    int next_surf2[]= {0, 2, 2, 0};

    int nq = 4; // number of quadrature currents. 
    int counter_j[] = {1, 0, 3, 2};

    for (int y = 0; y < ch; y++)
    {
        for (int x = 0; x < cw; x++)
        {
            int ii = y * cw + x;
            meshCell = &_cells[ii];

            for (int i = 0; i < 4; i++)
            {
                surfaceCorner = meshCell->getMeshSurfaces(i + 4);

                /* perform splitting inside of this cell: distribute
                 * evenly to the two surfaces surrounding the
                 * corner */
                surface1 = meshCell->getMeshSurfaces(surf[i]);
                surface2 = meshCell->getMeshSurfaces(surf2[i]);

                for (int j = 0; j < nq; j++)
                {
                    for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
                    {
                        current = f * surfaceCorner->getQuadCurrent(g, j);
                        surface1->incrementQuadCurrent(current, g, j);
                        surface2->incrementQuadCurrent(current, g, j);
                    }
                }

                /* effects on the neighboring cells */
                /* if left or right cell exists */
                if ((x > min_x[i]) && (x < max_x[i]))
                {
                    meshCellNext = 
                        &_cells[ii + next_x[i]];
                    surface_next1 = meshCellNext->getMeshSurfaces(next_surf[i]);
                   
                     for (int j = 0; j < nq; j++)
                     {
                         for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
                         {
                             current = f * surfaceCorner->getQuadCurrent(g, j);
                             surface_next1->incrementQuadCurrent(current, g, j);
                         }
                     }
                } 
                /* only need to worry about this when reflective */
                else if (((x == min_x[i]) && (_boundary[0] == REFLECTIVE)) 
                         || ((x == max_x[i]) && (_boundary[2] == REFLECTIVE)))
                {
                    surface_next1 = meshCell->getMeshSurfaces(next_surf[i]);
                   
                     for (int j = 0; j < nq; j++)
                     {
                         for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
                         {
                             current = f * surfaceCorner->getQuadCurrent(g, j);
                             surface_next1->incrementQuadCurrent
                                 (current, g, counter_j[j]);
                         }
                     }
                }



                if ((y > min_y[i]) && (y < max_y[i])) 
                {
                     meshCellNext2 = 
                        &_cells[ii + next_y2[i] * cw];
                     surface_next2 = meshCellNext2->getMeshSurfaces
                        (next_surf2[i]);

                     for (int j = 0; j < nq; j++)
                     {
                         for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
                         {
                             current = f * surfaceCorner->getQuadCurrent(g, j);
                             surface_next2->incrementQuadCurrent(current, g, j);
                         }
                     }
                } 
                else if (((y == min_y[i]) && (_boundary[3] == REFLECTIVE)) 
                         || ((y == max_y[i]) && (_boundary[1] == REFLECTIVE)))
                {
                     surface_next2 = meshCell->getMeshSurfaces(next_surf2[i]);

                     for (int j = 0; j < nq; j++)
                     {
                         for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
                         {
                             current = f * surfaceCorner->getQuadCurrent(g, j);
                             surface_next2->incrementQuadCurrent
                                 (current, g, counter_j[j]);
                         }
                     }
                } 



            }
        }
    }
    return;
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

void Mesh::computeTotQuadCurrents(){

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

            /* loop over index, energy groups */
            for (int j = 0; j < 2; j++)
            {
                /* set current tally to 0 */
                sum_cur = 0.0;

                for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
                    sum_cur += surfaceSide->getQuadCurrent(e, j);

                /* set current at group 0 to total current */
                surfaceSide->setQuadCurrent(sum_cur, 0, j);
            }
        }
    }
}

void Mesh::setKeffCMFD(double keff, int iter){
    _keff_cmfd[iter] = keff;
}

double Mesh::getKeffCMFD(int iter){
    return _keff_cmfd[iter];
}

void Mesh::setKeffMOC(double keff, int iter){
    _keff_moc[iter] = keff;
}

double Mesh::getKeffMOC(int iter){
    return _keff_moc[iter];
}

MeshSurface **Mesh::getSurfaces(){
    return _surfaces;
}

double Mesh::getOldTime(){
    return _old_time;
}

void Mesh::setOldTime(double time){
    _old_time = time;
}

/**
 * @brief Gets the mesh level
 * @return _mesh_level mesh level
 **/
int Mesh::getMeshLevel(){
    return _mesh_level;
}


/**
 * @brief Sets the cmfd level
 * @parap cmfd_level cmfd level
 **/
void Mesh::setMeshLevel(int mesh_level){
    _mesh_level = mesh_level;
}
