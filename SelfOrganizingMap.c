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

#include "SelfOrganizingMap.hxx"

#ifdef __cplusplus
extern "C" {
#endif

// method:
//	som_allocate
//
// description:
//
SelfOrganizingMap_t * som_allocate( int rows, int columns, int dimension )
{
	SelfOrganizingMap_t * som = (SelfOrganizingMap_t*)( malloc( sizeof( SelfOrganizingMap_t ) ) );

	// som attributes
	som->_rows = rows;
	som->_columns = columns;
	som->_dimension = dimension;

	// som state
	som->_decay = 0.0;
	som->_train_iterations = 0;
	som->_bmu_column = 0;
	som->_bmu_row = 0;
	som->_bmu_distance = 0.0;

	// allocate label and unit memory
	int size = rows * columns;
	som->_labels = (char**)( malloc( sizeof( char* ) * size ) );
	som->_units = (Vector_t**)( malloc( sizeof( Vector_t* ) * size ) );

	// initialize unit vectors
	double * data = (double*)( malloc( sizeof( double ) * dimension ) );
	for( int p = 0; p < size; p++ ) 
	{
		// generate a random gaussian vector
		for( int d = 0; d < dimension; d++ )
		{
			// generate a single random number from -1.0 .. 1.0
			data[ d ] = ( (double)( rand() ) / (double)( RAND_MAX ) ) * 2.0 - 1.0;
		}
		Vector_t * element = mat_allocate1d_with_data( dimension, data );

		// normalize to unit hypersphere
		mat_normalize( element );

		// store
		som->_units[ p ] = element;
		som->_labels[ p ] = "-";
	}

	// clean-up
	free( data );

	// return result
	return som;
}

// method:
//	som_deallocate
//
// description:
//
void som_deallocate( SelfOrganizingMap_t * som )
{
	assert( som != 0 );
	
	if( som != 0 )
	{
		int size = som->_rows * som->_columns;
		for( int p = 0; p < size; p++ )
		{
			mat_deallocate( som->_units[ p ] );
		}
		free( som->_units );
		free( som->_labels );
		free( som );
	}
}

// method:
//	som_train
//
// description:
//	Wv(t + 1) = Wv(t) + theta( v, t )alpha( t )( D(t) - Wv(t) )
//
void som_train( SelfOrganizingMap_t * som, Vector_t * observations[], int count, int iterations, double tolerance )
{
	som_train_with_labels( som, observations, NULL, count, iterations, tolerance );
}

// method:
//	som_train_with_labels
//
// description:
//	Wv(t + 1) = Wv(t) + theta( v, t )alpha( t )( D(t) - Wv(t) )
//
void som_train_with_labels( SelfOrganizingMap_t * som, Vector_t * observations[], char ** labels, int count, int iterations, double tolerance )
{
	assert( som != 0 );
	if( som == 0 ) return;

	// initialize trained iteration count
	som->_train_iterations = iterations;
	if( iterations == 0 )
	{
		som->_train_iterations = 1;
		if( ( tolerance < 1.0 ) && ( tolerance > 0.0 ) )
		{
			som->_train_iterations = 1.0 / tolerance;
		}
	}

	// train for the given number of iterations or tolerance
	double running_average_distance = 0.0;
	double old_running_average_distance = 0.0;
	for( int i = 0; i < som->_train_iterations; i++ )
	{
		for( int x = 0; x < count; x++ )
		{
			Vector_t * observation = observations[ x ];

			// find best matching unit and euclidian distance
			Vector_t * bmu = som_evaluate( som, observation );
			if( labels != NULL )
			{
				int position = som->_bmu_row * som->_columns + som->_bmu_column;
				som->_labels[ position ] = labels[ x ];
			}

			// update all nodes according to relative location and iteration
			double alpha = som_alpha( som, i );
			for( int r = 0; r < som->_rows; r++ )
			{
				for( int c = 0; c < som->_columns; c++ )
				{
					// create copy of observation for this iteration
					Vector_t * delta = mat_allocate_with_mat( observation );

					// multiply by neigborhood and learning coefficient
					double magnitude = som_theta( som, r, c, i ) * alpha;
					mat_subtract_mat( delta, bmu );
					mat_multiply( delta, magnitude );

					// set new value of bmu
					mat_add_mat( bmu, delta );

					// clean-up
					mat_deallocate( delta );
				}
			}

			// determine running average of distance
			double distance = som->_bmu_distance;
			double iteration = (double)(i + 1);
			running_average_distance = ( (running_average_distance * (iteration - 1.0)) + distance ) / iteration;
		}

		// stop training if the old running average isn't changeing more than the given tolerance
		if( abs( running_average_distance - old_running_average_distance ) < tolerance )
		{
			break;
		}
		old_running_average_distance = running_average_distance;
	}
}


// method:
//	som_evaluate
//
// description:
//
Vector_t * som_evaluate( SelfOrganizingMap_t * som, Vector_t * observation )
{
	assert( som != 0 );
	assert( observation != 0 );

	Vector_t * result = 0;
	double smallest_distance = -1.0;
	for( int r = 0; r < som->_rows; r++ )
	{
		for( int c = 0; c < som->_columns; c++ )
		{
			int p = r * som->_columns + c;
			Vector_t * v = som->_units[ p ];
			double distance = mat_distance( v, observation );
			if( ( smallest_distance == -1.0 ) || ( distance < smallest_distance ) )
			{
				smallest_distance = distance;
				som->_bmu_distance = distance;
				som->_bmu_row = r;
				som->_bmu_column = c;
				result = v;
			}
		}
	}
	return result;
}

// method:
//	som_evaluated_label
//
// description:
//
char * som_evaluated_label( SelfOrganizingMap_t * som )
{
	assert( som != 0 );

	char * result = 0;
	if( som != 0 )
	{
		int position = som->_bmu_row * som->_columns + som->_bmu_column;
		result = som->_labels[ position ];		
	}
	return result;
}

// method:
//	som_evaluated_bmu
//
// description:
//
Vector_t * som_evaluated_bmu( SelfOrganizingMap_t * som )
{
	assert( som != 0 );

	Vector_t * result = 0;
	if( som != 0 )
	{
		int position = som->_bmu_row * som->_columns + som->_bmu_column;
		result = som->_units[ position ];		
	}
	return result;
}

// method:
//	som_evaluated_bmu_position
//
// description:
//
void som_evaluated_bmu_position( SelfOrganizingMap_t * som, int * bmu_row, int * bmu_column )
{
	assert( som != 0 );
	assert( bmu_row != 0 );
	assert( bmu_column != 0 );
	
	if( ( som != 0 ) && ( bmu_row != 0 ) && ( bmu_column != 0 ) )
	{
		*(bmu_row) = som->_bmu_row;
		*(bmu_column) = som->_bmu_column;
	}
}

// method:
//	som_umatrix
//
// description:
//
Matrix_t * som_umatrix( SelfOrganizingMap_t * som )
{
	assert( som != 0 );
	
	Matrix_t * result = 0;
	if( som != 0 )
	{
		int umat_rows = som->_rows * 2;
		int umat_cols = som->_columns;
		result = mat_allocate2d( umat_rows, umat_cols );
		for( int r = 0; r < umat_rows; r++ )
		{
			for( int c = 0; c < umat_cols; c++ )
			{
				int phase = r % 2;
				int som_row_a = ( r - phase ) / 2;
				int som_col_a = c;
				int som_position_a = som_row_a * som->_columns + som_col_a;
				
				double distance = 0.0;
				if( phase == 0 )
				{
					// every even row, calculate the distance between columns
					int som_row_b = som_row_a;
					int som_col_b = (som_col_a + 1) % som->_columns;
					int som_position_b = som_row_b * som->_columns + som_col_b;
					
					Vector_t * a = som->_units[ som_position_a ];
					Vector_t * b = som->_units[ som_position_b ];
					distance = mat_distance( a, b );
				}
				else
				{
					// every odd row, calculate the distance between rows
					int som_row_b = (som_row_a + 1) % som->_rows;
					int som_col_b = som_col_a;
					int som_position_b = som_row_b * som->_columns + som_col_b;
					
					Vector_t * a = som->_units[ som_position_a ];
					Vector_t * b = som->_units[ som_position_b ];
					distance = mat_distance( a, b );					
				}
				mat_set2d( result, distance, r, c );
			}
		}
		
	}
	return result;
}
    
// method:
//	som_umatrix3d
//
// description:
//
Matrix_t * som_umatrix3d( SelfOrganizingMap_t * som )
{
    assert( som != 0 );
    
    Matrix_t * result = 0;
    if( som != 0 )
    {
        int umat_rows = som->_rows;
        int umat_cols = som->_columns;
        int umat_depth = 2;
        result = mat_allocate3d( umat_rows, umat_cols, umat_depth );
        for( int r = 0; r < umat_rows; r++ )
        {
            for( int c = 0; c < umat_cols; c++ )
            {
                int som_position_a = r * som->_columns + c;
                
                double distance = 0.0;
                
                // calculate the distance between columns
                {
                    int som_row_b = r;
                    int som_col_b = (c + 1) % som->_columns;
                    int som_position_b = som_row_b * som->_columns + som_col_b;
                    
                    Vector_t * a = som->_units[ som_position_a ];
                    Vector_t * b = som->_units[ som_position_b ];
                    distance = mat_distance( a, b );
                    
                    mat_set3d( result, distance, r, c, 0 );
                }
                // every odd row, calculate the distance between rows
                {
                    int som_row_b = (r + 1) % som->_rows;
                    int som_col_b = c;
                    int som_position_b = som_row_b * som->_columns + som_col_b;
                    
                    Vector_t * a = som->_units[ som_position_a ];
                    Vector_t * b = som->_units[ som_position_b ];
                    distance = mat_distance( a, b );					
                    
                    mat_set3d( result, distance, r, c, 1 );
                }
            }
        }
    }
    return result;
}

// method:
//	som_theta
//
// description:
//	neighborhood function
//	i.e., ( 1 / ::sqrt( 2 * pi ) ) * ::exp( ( -1 ) * ( 1 / 2 ) * ::pow( distance, 2 ) )
//  a.k.a., normal distrubution function
//
double som_theta( SelfOrganizingMap_t * som, int row, int column, int t )
{
	assert( som != 0 );

	double result = 0.0;
	if( som != 0 )
	{
		double row_diff_squared = (som->_bmu_row - row) * (som->_bmu_row - row);
		double col_diff_squared = (som->_bmu_column - column ) * (som->_bmu_column - column );
		double distance = sqrt( row_diff_squared + col_diff_squared );
		result = 0.39894228 * exp( -0.5 * distance * distance );
	}
	return result;
}

// method:
//	som_alpha
//
// description:
//	monotonically decreasing learning coefficient
//
double som_alpha( SelfOrganizingMap_t * som, int t )
{
	assert( som != 0 );
	assert( som->_train_iterations != 0.0 );

	double result = som->_decay;
    if( som->_decay == 0.0 )
    {
        if( ( som != 0 ) && ( som->_train_iterations != 0.0 ) )
        {
            result = ( som->_train_iterations - t ) / som->_train_iterations;
        }
    }
	return result;
}

#ifdef __cplusplus
}
#endif
