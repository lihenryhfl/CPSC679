//////////////////////////////////////////////////////////////////////
// This file is part of Wavelet Turbulence.
//
// Wavelet Turbulence is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Wavelet Turbulence is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Wavelet Turbulence.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2008 Theodore Kim and Nils Thuerey
//
//////////////////////////////////////////////////////////////////////
//
#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>

//////////////////////////////////////////////////////////////////////
// NT helper functions
//////////////////////////////////////////////////////////////////////
template < class T > inline T ABS( T a ) {
	return (0 < a) ? a : -a ;
}

template < class T > inline void SWAP_POINTERS( T &a, T &b ) {
	T temp = a;
	a = b;
	b = temp;
}

template < class T > inline void CLAMP( T &a, T b=0., T c=1.) {
	if(a<b) { a=b; return; }
	if(a>c) { a=c; return; }
}

template < class T > inline T MIN( T a, T b) {
	return (a < b) ? a : b;
}

template < class T > inline T MAX( T a, T b) {
	return (a > b) ? a : b;
}

template < class T > inline T MAX3( T a, T b, T c) {
	T max = (a > b) ? a : b;
	max = (max > c) ? max : c;
	return max;
}

template < class T > inline float MAX3V( T vec) {
	float max = (vec[0] > vec[1]) ? vec[0] : vec[1];
	max = (max > vec[2]) ? max : vec[2];
	return max;
}

template < class T > inline float MIN3V( T vec) {
	float min = (vec[0] < vec[1]) ? vec[0] : vec[1];
	min = (min < vec[2]) ? min : vec[2];
	return min;
}

#endif
