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
#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <iostream>

namespace INTERPOLATE {

//////////////////////////////////////////////////////////////////////
// linear interpolators
//////////////////////////////////////////////////////////////////////
static inline float lerp(float t, float a, float b) {
	return ( a + t * (b - a) );
}

static inline float lerp(float* field, float x, float y, int res) {
	// clamp backtrace to grid boundaries
	if (x < 0.5f) x = 0.5f;
	if (x > res - 1.5f) x = res - 1.5f;
	if (y < 0.5f) y = 0.5f;
	if (y > res - 1.5f) y = res - 1.5f;

	const int x0 = (int)x;
	const int y0 = (int)y;
	x -= x0;
	y -= y0;
	float d00, d10, d01, d11;

	// lerp the velocities
	d00 = field[x0 + y0 * res];
	d10 = field[(x0 + 1) + y0 * res];
	d01 = field[x0 + (y0 + 1) * res];
	d11 = field[(x0 + 1) + (y0 + 1) * res];
	return lerp(y, lerp(x, d00, d10),
			lerp(x, d01, d11));
}

//////////////////////////////////////////////////////////////////////////////////////////
// interpolate a vector from 2 fields, in 2D
//////////////////////////////////////////////////////////////////////////////////////////
static inline void lerp2dVec(FIELD_2D& field1, FIELD_2D& field2, float x, float y, float* out) {
  int xres = field1.xRes();
  int yres = field1.yRes();
  assert (xres == field2.xRes() && yres == field2.yRes());
  // clamp pos to grid boundaries
  if (x < 0.5) x = 0.5;
  if (x > xres - 1.5) x = xres - 1.5;
  if (y < 0.5) y = 0.5;
  if (y > yres - 1.5) y = yres - 1.5;

  // locate neighbors to interpolate
  const int x0 = (int) x;
  const int x1 = x0 + 1;
  const int y0 = (int) y;
  const int y1 = y0 + 1;

  // get interpolation weights
  const float s1 = x - (float)x0;
  const float s0 = 1.0f - s1;
  const float t1 = y - (float)y0;
  const float t0 = 1.0f - t1;

  out[0] = s0 * (t0 * field1(x0, y0) + t1 * field1(x0, y1)) +
    s1 * (t0 * field1(x1, y0) + t1 * field1(x1, y1));

  out[1] = s0 * (t0 * field2(x0, y0) + t1 * field2(x0, y1)) +
    s1 * (t0 * field2(x1, y0) + t1 * field2(x1, y1));
}

//////////////////////////////////////////////////////////////////////////////////////////
// interpolate a vector from 2 fields, in 2D
//////////////////////////////////////////////////////////////////////////////////////////
static inline float lerp2d(FIELD_2D& field, float x, float y) {
  int xres = field.xRes();
  int yres = field.yRes();

  // clamp pos to grid boundaries
  if (x < 0.5) x = 0.5;
  if (x > xres - 1.5) x = xres - 1.5;
  if (y < 0.5) y = 0.5;
  if (y > yres - 1.5) y = yres - 1.5;

  // locate neighbors to interpolate
  const int x0 = (int) x;
  const int x1 = x0 + 1;
  const int y0 = (int) y;
  const int y1 = y0 + 1;

  // get interpolation weights
  const float s1 = x - (float)x0;
  const float s0 = 1.0f - s1;
  const float t1 = y - (float)y0;
  const float t0 = 1.0f - t1;

  return s0 * (t0 * field(x0, y0) + t1 * field(x0, y1)) +
    s1 * (t0 * field(x1, y0) + t1 * field(x1, y1));
}

};
#endif
