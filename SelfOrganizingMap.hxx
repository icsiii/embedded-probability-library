#pragma once
//
// This file is part of the skunk-works/math package.
// 
// The skunk-works/math package is free software: you can redistribute it 
// and/or modify it under the terms of the GNU Lesser General Public 
// License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version.
// 
// The skunk-works/math package is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// See <http://www.gnu.org/licenses/> for more information.
//
// Author: Ira Stevens
//

#include "Matrix.hxx"

typedef struct
{
	int _train_iterations;

	// last calculated bmu position and distance
	double _bmu_distance;
	int _bmu_row;
	int _bmu_column;

	// monotonically decreasing learning coefficient decay
	double _decay; 

	// the self organizing map
	int _rows;				// unit and label size in rows
	int _columns;			// unit and label size in columns
	int _dimension;			// vector size of units
	Vector_t ** _units;		// unit vectors
	char ** _labels;		// unit labels
	void ** _data;			// unit data
} SelfOrganizingMap_t;

#ifdef __cplusplus
extern "C" {
#endif

SelfOrganizingMap_t * som_allocate( int rows, int columns, int dimension );
void		som_deallocate( SelfOrganizingMap_t * som );

void		som_train( SelfOrganizingMap_t * som, Vector_t * observations[], int count, int iterations, double tolerance );
void		som_train_with_labels( SelfOrganizingMap_t * som, Vector_t * observations[], char ** labels, int count, int iterations, double tolerance );
Vector_t *	som_evaluate( SelfOrganizingMap_t * som, Vector_t * observation );
char *		som_evaluated_label( SelfOrganizingMap_t * som );
Vector_t *	som_evaluated_bmu( SelfOrganizingMap_t * som );
void		som_evaluated_bmu_position( SelfOrganizingMap_t * som, int * bmu_row, int * bmu_column );
Matrix_t *	som_umatrix( SelfOrganizingMap_t * som );
Matrix_t *  som_umatrix3d( SelfOrganizingMap_t * som );

double		som_theta( SelfOrganizingMap_t * som, int row, int column, int t ); // neighborhood function
double		som_alpha( SelfOrganizingMap_t * som, int t ); // monotonically decreasing learning coefficient

#ifdef __cplusplus
}
#endif
