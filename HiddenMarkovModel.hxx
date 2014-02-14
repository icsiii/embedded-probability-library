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
// Author: This is a C port of LGPL HMM .NET implementation from 
//         Accord-NET, see github.com/accord-net/framework
//

#include "Matrix.hxx"

typedef enum 
{
	HiddenMarkovModelTypeErgodic,	// forward-backward
	HiddenMarkovModelTypeForward	// forward-only
} HiddenMarkovModelType_t;

typedef struct
{
	HiddenMarkovModelType_t _type;	// HMM type (forward-only or forward-backward)
	int _symbols;					// number of symbols per state
	int _states;					// number of state transistions
	Matrix_t * _pi;					// vector of initial states
	Matrix_t * _A;					// matrix (2d) of transisions
	Matrix_t * _B;					// matrix (2d) of emissions
	char * _label;					// associated label
	void * _data;					// associated data
} HiddenMarkovModel_t;

#ifdef __cplusplus
extern "C" {
#endif

HiddenMarkovModel_t * hmm_allocate( int symbols, int states, HiddenMarkovModelType_t type );
void		hmm_deallocate( HiddenMarkovModel_t * );

double		hmm_train( HiddenMarkovModel_t * hmm, Vector_t * observations[], int count, int iterations, double tolerance );
double		hmm_evaluate( HiddenMarkovModel_t * hmm, Vector_t * observation );

#ifdef __cplusplus
}
#endif
