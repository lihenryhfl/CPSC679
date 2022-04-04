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

// for e.g. MSVC compiler...
// some of these defines can be needed
// for linux systems as well (e.g. FLT_MAX)
#ifndef __FLT_MAX__
#	ifdef FLT_MAX  // try to use it instead
#		define __FLT_MAX__ FLT_MAX
#	else // FLT_MAX
#		define __FLT_MAX__ 3.402823466e+38f
#	endif // FLT_MAX
#endif // __FLT_MAX__
#ifndef __DBL_MAX__
#	ifdef DBL_MAX // try to use it instead
#		define __DBL_MAX__ DBL_MAX
#	else // DBL_MAX
#		define __DBL_MAX__ 1.7976931348623158e+308
#	endif // DBL_MAX
#endif // __DBL_MAX__

#ifndef FLT_MAX
#define FLT_MAX __FLT_MAX__
#endif

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
