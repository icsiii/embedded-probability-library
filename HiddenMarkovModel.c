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
// Author: This is a C port of LGPL HMM .NET implementation from 
//         Accord-NET, see github.com/accord-net/framework.
//

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "HiddenMarkovModel.hxx"

#ifdef __cplusplus
extern "C" {
#endif
    
// method:
//	hmm_allocate
//
// description:
//
HiddenMarkovModel_t * hmm_allocate( int symbols, int states, HiddenMarkovModelType_t type )
{
	assert( symbols > 0 );
	assert( states > 0 );
	
	HiddenMarkovModel_t * hmm = (HiddenMarkovModel_t*)( malloc( sizeof( HiddenMarkovModel_t ) ) );
	hmm->_symbols = symbols;
	hmm->_states = states;
	hmm->_type = type;
	
	// initialize initial state
	hmm->_pi = mat_allocate1d( states );
	mat_set1d( hmm->_pi, 1.0, 0 );
	
	// initialize transistions
	hmm->_A = mat_allocate2d( states, states );
	for( int i = 0; i < states; i++ )
	{
		for( int j = 0; j < states; j++ )
		{
			if( type == HiddenMarkovModelTypeErgodic )
			{
				mat_set2d( hmm->_A, (1.0 / (double)states), i, j );
			}
			else
			{
				mat_set2d( hmm->_A, (1.0 / (double)(states - i)), i, j );
			}
		}
	}
	
	// initialize emissions
	hmm->_B = mat_allocate2d( states, symbols );
	for( int i = 0; i < states; i++ )
	{
		for( int j = 0; j < symbols; j++ )
		{
			mat_set2d( hmm->_B, (1.0 / (double)symbols), i, j );
		}
	}
	
	// return result
	return hmm;
}

// method:
//	hmm_deallocate
//
// description:
//
void hmm_deallocate( HiddenMarkovModel_t * hmm )
{
	assert( hmm != 0 );
	if( hmm != 0 )
	{
		mat_deallocate( hmm->_B );
		mat_deallocate( hmm->_A );
		mat_deallocate( hmm->_pi );
		free( hmm );
	}
}

// method:
//	hmm_forward( hmm : HMM*, sequence : Vector*, scaling : Vector** ) : Matrix*
//
// description:
//	the baum-welch forward pass with scaling
//
// return:
//	the returned parameter points to memory owned by the user, it should be 
//	deallocted using mat_deallocate( ... ).
//
// parameter hmm:
//	the object of the method call
//
// parameter sequence:
//	the observation sequence
//
// parameter scaling:
//	the returned scaling vector, the caller must deallocate the memory
//	returned using mat_deallocate( ... ).
//
static Matrix_t * hmm_forward( HiddenMarkovModel_t * hmm, Vector_t * sequence, Vector_t ** scaling )
{
	assert( hmm != 0 );
	int T = mat_xsize( sequence );
	Matrix_t * fwd = mat_allocate2d( T, hmm->_states );
	
	// ** NOTE ** scaling vector is initialized to all zeros at this point
	// ** NOTE ** scaling vector overwritten with a new pointer (user beware)
	*scaling = mat_allocate1d( T );
	Vector_t * c = *scaling;
	
	// initialization
	for( int i = 0; i < hmm->_states; i++ )
	{
		int ov = (int)(mat_get1d( sequence, 0 ));
		double fwdv = mat_get1d( hmm->_pi, i ) * mat_get2d( hmm->_B, i, ov );
		double cv = mat_get1d( c, 0 );
		mat_set2d( fwd, fwdv, 0, i );
		mat_set1d( c, cv + fwdv, 0 );
	}
	double cv = mat_get1d( c, 0 );
	if( cv != 0.0 )
	{
		for( int i = 0; i < hmm->_states; i++ )
		{
			double fwdv = mat_get2d( fwd, 0, i ) / mat_get1d( c, 0 );
			mat_set2d( fwd, fwdv, 0, i );
		}
	}
	
	// inducation
	for( int t = 1; t < T; t++ )
	{
		for( int i = 0; i < hmm->_states; i++ )
		{
			double ov = (int)(mat_get1d( sequence, t ));
			double p = mat_get2d( hmm->_B, i, ov );
			
			double sum = 0.0;
			for( int j = 0; j < hmm->_states; j++ )
			{
				sum += mat_get2d( fwd, t - 1, j ) * mat_get2d( hmm->_A, j, i );
			}
			
			double fwdv = sum * p;
			mat_set2d( fwd, fwdv, t, i );

			double cv = mat_get1d( c, t );
			mat_set1d( c, cv + fwdv, t );
		}
		double cv = mat_get1d( c, t );
		if( cv != 0 )
		{
			for( int i = 0; i < hmm->_states; i++ )
			{
				double fwdv = mat_get2d( fwd, t, i ) / mat_get1d( c, t );
				mat_set2d( fwd, fwdv, t, i );
			}
		}
	}
	
	// return result
	return fwd;
}

// method:
//	hmm_backward( hmm : HMM*, sequence : Vector*, scaling : Vector** ) : Matrix*
//
// description:
//	the baum-welch backward pass with scaling
//
// return:
//	the returned parameter points to memory owned by the user, it should be 
//	deallocted using mat_deallocate( ... ).
//
// parameter hmm:
//	the object of the method call
//
// parameter sequence:
//	the observation sequence
//
// parameter scaling:
//	the scaling vector created from the forward method
//
static Matrix_t * hmm_backward( HiddenMarkovModel_t * hmm, Vector_t * sequence, Vector_t * scaling )
{
	assert( hmm != 0 );
	int T = mat_xsize( sequence );
	Matrix_t * bwd = mat_allocate2d( T, hmm->_states );
	
	Vector_t * c = scaling;
	
	// initialize
	for( int i = 0; i < hmm->_states; i++ )
	{
		double cv = mat_get1d( c, T - 1 );
		double bwdv = 1.0 / cv;
		mat_set2d( bwd, bwdv, T - 1, i );
	}
	
	// induction
	for( int t = T - 2; t >= 0; t-- )
	{
		for( int i = 0; i < hmm->_states; i++ )
		{
			double sum = 0.0;
			for( int j = 0; j < hmm->_states; j++ )
			{
				int ov = (int)(mat_get1d( sequence, t + 1 ));
				sum += mat_get2d( hmm->_A, i, j ) * mat_get2d( hmm->_B, j, ov ) * mat_get2d( bwd, t + 1, j );
			}
			double bwdv = mat_get2d( bwd, t, i ) + sum / mat_get1d( c, t );
			mat_set2d( bwd, bwdv, t, i );
		}
	}
	
	// return result
	return bwd;
}

// method:
//	hmm_has_converged( hmm : HMM*, old_likelihood : double, new_likelihood : double, current_iteration : int, max_iteration : int, tolerance : double ) : int
//
// description:
//	indicate if the hmm has converged within the given tolerance
//
static int hmm_has_converged( double old_likelihood, double new_likelihood, int current_iteration, int max_iteration, double tolerance )
{
	int result = 0;
	if( isinf( new_likelihood ) || isnan( new_likelihood ) )
	{
		result = 1;
	}
	if( ( tolerance > 0.0 ) && ( fabs( old_likelihood - new_likelihood ) < tolerance ) )
	{
		result = 1;
	}
	if( ( max_iteration > 0 ) && ( current_iteration >= max_iteration ) )
	{
		result = 1;
	}
	return result;
}

// method:
//	hmm_train
//
// description:
//
double hmm_train( HiddenMarkovModel_t * hmm, Vector_t * observations[], int count, int iterations, double tolerance )
{
	assert( hmm != 0 );
	double new_likelihood = 0.0;
	if( ( iterations != 0 ) || ( tolerance != 0.0 ) )
	{
		//int N = sizeof( observations ) / sizeof( Vector_t * );
		int N = count;
		int current_iteration = 1;
		int stop = 0;
		
		// initialize epsilon (aka, ksi or psi) and gamma
		Matrix_t * epsilon[ N ];
		Matrix_t * gamma[ N ];
		for( int i = 0; i < N; i++ )
		{
			int T = mat_xsize( observations[ i ] );
			epsilon[ i ] = mat_allocate3d( T, hmm->_states, hmm->_states );
			gamma[ i ] = mat_allocate2d( T, hmm->_states );
		}

		// initial log likelihood
		double old_likelihood = 0.0;
		
		// train until done (max iterations or converged within tolerance)
		do
		{
			// train for each sequence in observations
			for( int i = 0; i < N; i++ )
			{
				Vector_t * sequence = observations[ i ];
				int T = mat_xsize( sequence );
				Vector_t * scaling = 0;
				
				// (a) calculate forward and backward probability
				Matrix_t * fwd = hmm_forward( hmm, sequence, &scaling );
				Matrix_t * bwd = hmm_backward( hmm, sequence, scaling );
				
				// (b) calculate the frequency of the transition-emission pair valus
				//     and divide by the probability of the entire sequence
				//
				// calculate gamma
				for( int t = 0; t < T; t++ )
				{
					double s = 0.0;
					for( int k = 0; k < hmm->_states; k++ )
					{
						double gv = mat_get2d( fwd, t, k ) * mat_get2d( bwd, t, k );
						mat_set2d( gamma[ i ], gv, t, k );
						s += gv;
					}
					if( s != 0.0 )
					{
						for( int k = 0; k < hmm->_states; k++ )
						{
							double gv = mat_get2d( gamma[ i ], t, k );
							mat_set2d( gamma[ i ], (gv / s), t, k );
						}
					}
				}
				
				// calculate epsilon 
				for( int t = 0; t < T - 1; t++ )
				{
					double s = 0.0;
					for( int k = 0; k < hmm->_states; k++ )
					{
						for( int l = 0; l < hmm->_states; l++ )
						{
							int next_symbol = (int)(mat_get1d( sequence, t + 1 ));
							double gv = mat_get2d( fwd, t, k ) * mat_get2d( bwd, t + 1, l );
							double ev = gv * mat_get2d( hmm->_A, k, l ) * mat_get2d( hmm->_B, l, next_symbol );
							mat_set3d( epsilon[ i ], ev, t, k, l );
							s += ev;
						}
					}
					if( s != 0.0 )
					{
						for( int k = 0; k < hmm->_states; k++ )
						{
							for( int l = 0; l < hmm->_states; l++ )
							{
								double ev = mat_get3d( epsilon[ i ], t, k, l );
								mat_set3d( epsilon[ i ], (ev / s ), t, k, l );
							}
						}
					}
				}
				
				// calculate log likelihood
				for( int t = 0; t < mat_xsize( scaling ); t++ )
				{
					new_likelihood += log( mat_get1d( scaling, t ) );
				}
				
				// free working fwd, bwd and scaling matrix
				mat_deallocate( fwd );
				mat_deallocate( bwd );
				mat_deallocate( scaling );
				scaling = 0;
			}
		
			// average likelihood
			new_likelihood /= (double)N;
			
			// check for convergence
			if( hmm_has_converged( old_likelihood, new_likelihood, current_iteration, iterations, tolerance ) != 0 )
			{
				stop = 1;
			}
			else
			{
				// (c) calculate parameter re-estimation
				++current_iteration;
				old_likelihood = new_likelihood;
				new_likelihood = 0.0;
				
				// re-estimate initial state
				for( int k = 0; k < hmm->_states; k++ )
				{
					double s = 0.0;
					for( int i = 0; i < N; i++ )
					{
						s += mat_get2d( gamma[ i ], 0, k );
					}
					mat_set1d( hmm->_pi, (s / N), k );
				}
				
				// re-estimate transition probabilities
				for( int i = 0; i < hmm->_states; i++ )
				{
					for( int j = 0; j < hmm->_states; j++ )
					{
						double den = 0.0;
						double num = 0.0;
						for( int k = 0; k < N; k++ )
						{
							int T = mat_xsize( observations[ k ] );
							for( int l = 0; l < T - 1; l++ )
							{
								double ev = mat_get3d( epsilon[ k ], l, i, j );
								double gv = mat_get2d( gamma[ k ], l, i );
								num += ev;
								den += gv;
							}
						}
						double av = (den != 0.0) ? num / den : 0.0;
						mat_set2d( hmm->_A, av, i, j );
					}
				}
				
				// re-estimation emission probabilities
				for( int i = 0; i < hmm->_states; i++ )
				{
					for( int j = 0; j < hmm->_symbols; j++ )
					{
						double den = 0.0;
						double num = 0.0;
						for( int k = 0; k < N; k++ )
						{
							int T = mat_xsize( observations[ k ] );
							for( int l = 0; l < T; l++ )
							{
								double gv = mat_get2d( gamma[ k ], l, i );
								int ov = (int)(mat_get1d( observations[ k ], l ));
								if( ov == j ) num += gv;
								den += gv;
							}
						}
						double bv = (num == 0.0) ? 1e-10 : num / den;
						mat_set2d( hmm->_B, bv, i, j );
					}
				}
			}
		}
		while( stop == 0 );
		
		// free epsilon and gamma
		for( int i = 0; i < N; i++ )
		{
			mat_deallocate( epsilon[ i ] );
			mat_deallocate( gamma[ i ] );
		}
	}
	return new_likelihood;
}

// method:
//	hmm_evaluate
//
// description:
//
double hmm_evaluate( HiddenMarkovModel_t * hmm, Vector_t * observation )
{
	assert( hmm != 0 );
	double likelihood = 0.0;
	Vector_t * coefficients = 0;
	
	// calculate forward probabilities
	Matrix_t * fwd = hmm_forward( hmm, observation, &coefficients );
	for( int i = 0; i < mat_xsize( coefficients ); i++ )
	{
		likelihood += log( mat_get1d( coefficients, i ) );
	}
	
	// clean up memory
	mat_deallocate( fwd );
	mat_deallocate( coefficients );
	
	// return result
	return exp( likelihood );
}
    
#ifdef __cplusplus
}
#endif
