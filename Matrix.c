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

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.hxx"

#ifdef __cplusplus
extern "C" {
#endif

// method:
//	mat_allocate3d
//
// description:
//
Matrix_t * mat_allocate3d( int xsize, int ysize, int zsize )
{
	assert( zsize > 0 );
	assert( ysize > 0 );
	assert( xsize > 0 );
	
	Matrix_t * m = (Matrix_t*)( malloc( sizeof( Matrix_t ) ) );
	m->_xsize = xsize;
	m->_ysize = ysize;
	m->_zsize = zsize;
	m->_magnitude = 1.0;

	m->_type = MatrixTypeMatrix3d;
	if( zsize == 1 )
	{
		if( ysize == 1 )
		{
			m->_type = MatrixTypeVector;
		}
		else
		{
			m->_type = MatrixTypeMatrix2d;
		}
	}
	
	m->_zdata = (double ***)( malloc( zsize * sizeof( double** ) ) );
	for( int z = 0; z < zsize; z++ )
	{
		m->_zdata[ z ] = (double **)( malloc( ysize * sizeof( double* ) ) );
		for( int y = 0; y < ysize; y++ )
		{
			m->_zdata[ z ][ y ] = (double *)(malloc( xsize * sizeof( double ) ) );
			for( int x = 0; x < xsize; x++ )
			{
				m->_zdata[ z ][ y ][ x ] = 0.0;
			}
		}
	}
	return m;
}

// method:
//	mat_allocate2d
//
// description:
//
Matrix_t * mat_allocate2d( int xsize, int ysize )
{
	return mat_allocate3d( xsize, ysize, 1 );
}

// method:
//	mat_allocate1d
//
// description:
//
Matrix_t * mat_allocate1d( int xsize )
{
	return mat_allocate3d( xsize, 1, 1 );
}

// method:
//	mat_allocate1d_with_data
//
// description:
//
Matrix_t * mat_allocate1d_with_data( int xsize, double data[] )
{
	Matrix_t * result = mat_allocate3d( xsize, 1, 1 );
	for( int x = 0; x < xsize; x++ )
	{
		mat_set1d( result, data[ x ], x );
	}
	return result;
}

// method:
//	mat_allocate_with_mat
//
// description:
//
Matrix_t * mat_allocate_with_mat( Matrix_t * v )
{
	assert( v != 0 );
	Matrix_t * result = mat_allocate3d( v->_xsize, v->_ysize, v->_zsize );
	mat_copy( result, v );
	return result;
}

// method:
//	mat_deallocate
//
// description:
//
void mat_deallocate( Matrix_t * m )
{
	assert( m != 0 );
	if( m != 0 )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				free( m->_zdata[ z ][ y ] );
			}
			free( m->_zdata[ z ] );
		}
		free( m->_zdata );
		free( m );
	}
}

// method:
//	mat_get3d
//
// description:
//
double mat_get3d( Matrix_t * m, int x, int y, int z )
{
	assert( m != 0 );
	double result = 0.0;
	if( ( x < m->_xsize ) && ( x >= 0 ) &&
		( y < m->_ysize ) && ( y >= 0 ) &&
		( z < m->_zsize ) && ( z >= 0 ) )
	{
		result = m->_zdata[ z ][ y ][ x ];
	}
	return result;
}

// method:
//	mat_get2d
//
// description:
//
double mat_get2d( Matrix_t * m, int x, int y )
{
	return mat_get3d( m, x, y, 0 );
}

// method:
//	mat_get1d
//
// description:
//
double mat_get1d( Matrix_t * m, int x )
{
	return mat_get3d( m, x, 0, 0 );
}

// method:
//	mat_set3d
//
// description:
//
void mat_set3d( Matrix_t * m, double v, int x, int y, int z )
{
	assert( m != 0 );
	if( ( x < m->_xsize ) && ( x >= 0 ) &&
	   ( y < m->_ysize ) && ( y >= 0 ) &&
	   ( z < m->_zsize ) && ( z >= 0 ) )
	{
		m->_zdata[ z ][ y ][ x ] = v;
	}
}

// method:
//	mat_set2d
//
// description:
//
void mat_set2d( Matrix_t * m, double v, int x, int y )
{
	mat_set3d( m, v, x, y, 0 );
}

// method:
//	mat_set1d
//
// description:
//
void mat_set1d( Matrix_t * m, double v, int x )
{
	mat_set3d( m, v, x, 0, 0 );
}

// method:
//	mat_xsize
//
// description:
//
int mat_xsize( Matrix_t * m )
{
	assert( m != 0 );
	return m->_xsize;
}

// method:
//	mat_ysize
//
// description:
//
int mat_ysize( Matrix_t * m )
{
	assert( m != 0 );
	return m->_ysize;
}

// method:
//	mat_zsize
//
// description:
//
int mat_zsize( Matrix_t * m )
{
	assert( m != 0 );
	return m->_zsize;
}

// method:
//	mat_copy
//
// description:
//
void mat_copy( Matrix_t * m, Matrix_t * v )
{
	assert( m != 0 );
	assert( v != 0 );
	assert( m != v );
	if( ( m != 0 ) && ( v != 0 ) && ( m != v ) )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = v->_zdata[ z ][ y ][ x ];
				}
			}
		}
	}
}

// method:
//	mat_maximum
//
// description:
//	calculate matrix maximum
//
double mat_maximum( Matrix_t * m )
{
	double result = 0.0;
	assert( m != 0 );
	if( m != 0 )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
                    if( m->_zdata[ z ][ y ][ x ] > result )
                    {
                        result = m->_zdata[ z ][ y ][ x ];
                    }
				}
			}
		}
	}
	return result;
}
    
// method:
//	mat_magnitude
//
// description:
//	calculate matrix magnitude (root-sum-square)
//
double mat_magnitude( Matrix_t * m )
{
	double result = 0.0;
	assert( m != 0 );
	if( m != 0 )
	{
		double sum = 0.0;
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					double square = m->_zdata[ z ][ y ][ x ] * m->_zdata[ z ][ y ][ x ];
					sum = sum + square;
				}
			}
		}
		result = sqrt( sum );
	}
	return result;
}

// method:
//	mat_normalize
//
// description:
//	normalize to unit vector
//
void mat_normalize( Matrix_t * m )
{
	m->_magnitude = mat_magnitude( m );
	mat_multiply( m, 1.0 / m->_magnitude );
}

// method:
//	mat_denormalize
//
// description:
//	denormalize the unit vector
//
void mat_denormalize( Matrix_t * m )
{
	mat_multiply( m, m->_magnitude );
	m->_magnitude = 1.0;
}

// method:
//	mat_distance
//
// description:
//	calculate euclidian distance
//
double mat_distance( Matrix_t * m, Matrix_t * v )
{
	double result = 0.0;
	assert( m != 0 );
	assert( m->_type == v->_type );
	assert( m->_xsize == v->_xsize );
	assert( m->_ysize == v->_ysize );
	assert( m->_zsize == v->_zsize );
	if( ( m != 0 ) && ( m->_type == v->_type ) && ( m->_xsize == v->_xsize ) && ( m->_ysize == v->_ysize ) && ( m->_zsize == v->_zsize ) )
	{
		double sum = 0.0;
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					double difference = m->_zdata[ z ][ y ][ x ] - v->_zdata[ z ][ y ][ x ];
					double difference_squared = difference * difference;
					sum = sum + difference_squared;
				}
			}
		}
		result = sqrt( sum );
	}
	return result;
}

// method:
//	mat_multiply
//
// description:
//
void mat_multiply( Matrix_t * m, double v )
{
	assert( m != 0 );
	if( m != 0 )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = m->_zdata[ z ][ y ][ x ] * v;
				}
			}
		}
	}
}

// method:
//	mat_transpose
//
// description:
//	a simple 1d/2d transpose. note that a 3d matrix will remain unchanged
//	in the z-axis (i.e., only the x and y axis are transposed.)
//
Matrix_t * mat_transpose( Matrix_t * m )
{
	assert( m != 0 );
	
	Matrix_t * result = 0;
	if( m != 0 ) 
	{
		result = mat_allocate3d( m->_ysize, m->_xsize, m->_zsize );
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					result->_zdata[ z ][ x ][ y ] = m->_zdata[ z ][ y ][ x ];
				}
			}
		}
	}
	return result;
}

// method:
//	mat_multiply_mat
//
// description:
//
void mat_multiply_mat( Matrix_t * m, Matrix_t * v )
{
	assert( m != 0 );
	assert( m->_type == v->_type );
	assert( m->_xsize == v->_xsize );
	assert( m->_ysize == v->_ysize );
	assert( m->_zsize == v->_zsize );
	if( ( m != 0 ) && ( m->_type == v->_type ) && ( m->_xsize == v->_xsize ) && ( m->_ysize == v->_ysize ) && ( m->_zsize == v->_zsize ) )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = m->_zdata[ z ][ y ][ x ] * v->_zdata[ z ][ y ][ x ];
				}
			}
		}
	}
}

// method:
//	mat_cross
//
// description:
//
void mat_cross( Matrix_t * m, Matrix_t * v )
{
	// ** incomplete **
}

// method:
//	mat_dot
//
// description:
//
void mat_dot( Matrix_t * m, Matrix_t * v )
{
	// ** incomplete **
}

// method:
//	mat_add
//
// description:
//
void mat_add( Matrix_t * m, double v )
{
	assert( m != 0 );
	if( m != 0 )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = m->_zdata[ z ][ y ][ x ] + v;
				}
			}
		}
	}
}

// method:
//	mat_add
//
// description:
//
void mat_add_mat( Matrix_t * m, Matrix_t * v )
{
	assert( m != 0 );
	assert( m->_type == v->_type );
	assert( m->_xsize == v->_xsize );
	assert( m->_ysize == v->_ysize );
	assert( m->_zsize == v->_zsize );
	if( ( m != 0 ) && ( m->_type == v->_type ) && ( m->_xsize == v->_xsize ) && ( m->_ysize == v->_ysize ) && ( m->_zsize == v->_zsize ) )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = m->_zdata[ z ][ y ][ x ] + v->_zdata[ z ][ y ][ x ];
				}
			}
		}
	}
}

// method:
//	mat_subtract
//
// description:
//
void mat_subtract( Matrix_t * m, double v )
{
	return mat_add( m, (-1.0) * v );
}

// method:
//	mat_subtract
//
// description:
//
void mat_subtract_mat( Matrix_t * m, Matrix_t * v )
{
	assert( m != 0 );
	assert( m->_type == v->_type );
	assert( m->_xsize == v->_xsize );
	assert( m->_ysize == v->_ysize );
	assert( m->_zsize == v->_zsize );
	if( ( m != 0 ) && ( m->_type == v->_type ) && ( m->_xsize == v->_xsize ) && ( m->_ysize == v->_ysize ) && ( m->_zsize == v->_zsize ) )
	{
		for( int z = 0; z < m->_zsize; z++ )
		{
			for( int y = 0; y < m->_ysize; y++ )
			{
				for( int x = 0; x < m->_xsize; x++ )
				{
					m->_zdata[ z ][ y ][ x ] = m->_zdata[ z ][ y ][ x ] - v->_zdata[ z ][ y ][ x ];
				}
			}
		}
	}
}

#ifdef __cplusplus
}
#endif
