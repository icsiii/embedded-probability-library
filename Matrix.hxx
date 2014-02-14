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

typedef enum
{
	MatrixTypeVector,
	MatrixTypeMatrix2d,
	MatrixTypeMatrix3d
} MatrixType_t;

typedef struct
{
	MatrixType_t _type;				// matrix type (3d, 2d or vector)
	double _magnitude;				// magnitude set when normalized, unset when denormalized
	int _xsize;						// x-dimension size
	int _ysize;						// y-dimension size
	int _zsize;						// z-dimension size
	double *** _zdata;				// the matrix elements
} Matrix_t;

typedef Matrix_t Vector_t;

#ifdef __cplusplus
extern "C" {
#endif

Matrix_t *	mat_allocate3d( int xsize, int ysize, int zsize );
Matrix_t *	mat_allocate2d( int xsize, int ysize );
Matrix_t *	mat_allocate1d( int xsize );
Matrix_t *	mat_allocate1d_with_data( int xsize, double data[] );
Matrix_t *	mat_allocate_with_mat( Matrix_t * v );
void		mat_deallocate( Matrix_t * );

double		mat_get3d( Matrix_t * m, int x, int y, int z );
double		mat_get2d( Matrix_t * m, int x, int y );
double		mat_get1d( Matrix_t * m, int x );

void		mat_set3d( Matrix_t * m, double v, int x, int y, int z );
void		mat_set2d( Matrix_t * m, double v, int x, int y );
void		mat_set1d( Matrix_t * m, double v, int x );

int			mat_xsize( Matrix_t * m );
int			mat_ysize( Matrix_t * m );
int			mat_zsize( Matrix_t * m );

void		mat_copy( Matrix_t * m, Matrix_t * v );

double      mat_maximum( Matrix_t * m );
double		mat_magnitude( Matrix_t * m );
void		mat_normalize( Matrix_t * m );
void		mat_denormalize( Matrix_t * m );
double		mat_distance( Matrix_t * m, Matrix_t * v );
void		mat_multiply( Matrix_t * m, double v );
void		mat_multiply_mat( Matrix_t * m, Matrix_t * v );
Matrix_t *	mat_transpose( Matrix_t * m );
//void		mat_cross( Matrix_t * m, Matrix_t * v );
//void		mat_dot( Matrix_t * m, Matrix_t * v );
void		mat_add( Matrix_t * m, double v );
void		mat_add_mat( Matrix_t * m, Matrix_t * v );
void		mat_subtract( Matrix_t * m, double v );
void		mat_subtract_mat( Matrix_t * m, Matrix_t * v );

#ifdef __cplusplus
}
#endif
