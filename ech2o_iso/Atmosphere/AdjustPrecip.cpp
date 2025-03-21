/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta
 *******************************************************************************/
/*
 * AdjustPrecip.cpp
 *
 *  Created on: Aug 20, 2010
 *      Author: Marco.Maneta
 */

#include "Atmosphere.h"

int Atmosphere::AdjustPrecip(float cdt){

	UINT4 r, c;


	for(UINT4 i = 0; i < _vSortedGrid.size(); i++)
		for (UINT4 j = 0; j < _vSortedGrid[i].cells.size() ; j++)
		{
			r = _vSortedGrid[i].cells[j].row;
			c = _vSortedGrid[i].cells[j].col;
			_PrecipDepth->matrix[r][c] = _Precip->matrix[r][c]* cdt; //yangx -for direct m output
			_Precip->matrix[r][c] *= _isohyet->matrix[r][c];

		}


	return EXIT_SUCCESS;

}
