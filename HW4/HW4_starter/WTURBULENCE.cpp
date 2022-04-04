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
// WTURBULENCE handling
///////////////////////////////////////////////////////////////////////////////////

#include "WTURBULENCE.h"
#include "INTERPOLATE.h"
#include "UTIL.h"
#include <MERSENNETWISTER.h>
#include "WAVELET_NOISE.h"
#include "EIGENVALUE_HELPER.h"
#include "LU_HELPER.h"
#include "FIELD_2D.h"

// needed to access static advection functions
#include "FLUID_2D_BOUNDED.h"

// 2^ {-5/6}
static const float persistence = 0.56123f;

//////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////
WTURBULENCE::WTURBULENCE(int xResSm, int yResSm, int amplify) :
  _densityBig(xResSm * amplify, yResSm * amplify),
  _densityBigOld(xResSm * amplify, yResSm * amplify),
  _bigUx(xResSm * amplify, yResSm * amplify),
  _bigUy(xResSm * amplify, yResSm * amplify),
  _tcU(xResSm, yResSm),
  _tcV(xResSm, yResSm),
  _tcTemp(xResSm, yResSm),
  _energy(xResSm, yResSm),
  _highFreqEnergy(xResSm, yResSm),
  _eigMin(xResSm, yResSm),
  _eigMax(xResSm, yResSm),
  _noiseTile(noiseTileSize, noiseTileSize),
  _resSm(xResSm, yResSm)
{
  // if noise magnitude is below this threshold, its contribution
  // is negilgible, so stop evaluating new octaves
  _cullingThreshold = 1e-3;

  // factor by which to increase the simulation resolution
  _amplify = amplify;

  // manually adjust the overall amount of turbulence
  _strength = 2.;

  // add the corresponding octaves of noise
  _octaves = log((float)_amplify) / log(2.0f);

  // noise resolution
  _xResBig = _amplify * xResSm;
  _yResBig = _amplify * yResSm;

  // original / small resolution
  _xResSm = xResSm;
  _yResSm = yResSm;

  // allocate high resolution density field
  _totalStepsBig = 0;

  // map all
  const float dx = 1./(float)(_resSm[0]);
  const float dy = 1./(float)(_resSm[1]);
  for (int y = 0; y < _yResSm; y++)
    for (int x = 0; x < _xResSm; x++)
    {
      _tcU(x, y) = x * dx;
      _tcV(x, y) = y * dy;
    }

  // noise tiles
  std::string noiseTileFilename = std::string("noise.wavelets");
  //std::string noiseTileFilename = std::string("noise.fft");
  generateTile(_noiseTile.data(), noiseTileFilename);
}

//////////////////////////////////////////////////////////////////////
// destructor
//////////////////////////////////////////////////////////////////////
WTURBULENCE::~WTURBULENCE() {
  delete[] _densityBig;
  delete[] _densityBigOld;

  delete[] _bigUx;
  delete[] _bigUy;

  delete[] _tcU;
  delete[] _tcV;
  delete[] _tcTemp;

  delete[] _eigMin;
  delete[] _eigMax;
  delete[] _noiseTile;

  delete[] _energy;
  delete[] _highFreqEnergy;
}

//////////////////////////////////////////////////////////////////////
// Get the smallest valid x derivative
//
// Takes the one-sided finite difference in both directions and
// selects the smaller of the two
//////////////////////////////////////////////////////////////////////
static float minDx(int x, int y, FIELD_2D& input)
{
  const int N = input.xRes() - 2;

  // get grid values
  float center = input(x, y);
  float left  = (x <= 1)    ? FLT_MAX : input(x - 1, y);
  float right = (x >= maxx) ? FLT_MAX : input(x + 1, y);

  const float dx = xRes;

  // get all the derivative estimates
  float dLeft   = (x <= 1)     ? FLT_MAX : (center - left) * dx;
  float dRight  = (x >= N)  ? FLT_MAX : (right - center) * dx;
  float dCenter = (x <= 1 || x >= N) ? FLT_MAX : (right - left) * dx * 0.5f;

  // if it's on a boundary, only one estimate is valid
  if (x <= 1) return dRight;
  if (x >= N) return dLeft;

  // if it's not on a boundary, get the smallest one
  float finalD;
  finalD = (fabs(dCenter) < fabs(dRight)) ? dCenter : dRight;
  finalD = (fabs(finalD)  < fabs(dLeft))  ? finalD  : dLeft;

  return finalD;
}

//////////////////////////////////////////////////////////////////////
// get the smallest valid y derivative
//
// Takes the one-sided finite difference in both directions and
// selects the smaller of the two
//////////////////////////////////////////////////////////////////////
static float minDy(int x, int y, FIELD_2D& input)
{
  const int N = input.yRes() - 2;

  // get grid values
  float center = input[index];
  float down  = (y <= 1) ? FLT_MAX : input(x, y - 1);
  float up = (y >= maxy) ? FLT_MAX : input(x, y + 1);

  const float dx = input.yRes(); // only for square domains

  // get all the derivative estimates
  float dDown   = (y <= 1)  ? FLT_MAX : (center - down) * dx;
  float dUp  = (y >= N)  ? FLT_MAX : (up - center) * dx;
  float dCenter = (y <= 1 || y >= N) ? FLT_MAX : (up - down) * dx * 0.5f;

  // if it's on a boundary, only one estimate is valid
  if (y <= 1) return dUp;
  if (y >= N) return dDown;

  // if it's not on a boundary, get the smallest one
  float finalD = (fabs(dCenter) < fabs(dUp)) ? dCenter : dUp;
  finalD = (fabs(finalD) < fabs(dDown)) ? finalD : dDown;

  return finalD;
}

//////////////////////////////////////////////////////////////////////
// handle texture coordinates (advection, reset, eigenvalues),
// Beware -- uses big density maccormack as temporary arrays
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::advectTextureCoordinates (float dtOrg, FIELD_2D& xvel, FIELD_2D& yvel) {
  // advection
  FLUID_2D_BOUNDED::swapFields(_tcTemp, _tcU);
  FLUID_2D_BOUNDED::copyBoundary(_tcTemp, true, true, 1);
  FLUID_2D_BOUNDED::advect(_tcU, _tcTemp, xvel, yvel);

  FLUID_2D_BOUNDED::swapFields(_tcTemp, _tcV);
  FLUID_2D_BOUNDED::copyBoundary(_tcTemp, true, true, 1);
  FLUID_2D_BOUNDED::advect(_tcV, _tcTemp, xvel, yvel);
}

//////////////////////////////////////////////////////////////////////
// Compute the eigenvalues of the advected texture
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::computeEigenvalues() {
  // stats
  float maxeig = -1.;
  float mineig = 10.;

  // texture coordinate eigenvalues
  for (int y = 1; y < _yResSm - 1; y++)
    for (int x = 1; x < _xResSm - 1; x++)
    {
      // compute jacobian
      float jacobian[2][2] = {
        { minDx(x, y, _tcU), minDx(x, y, _tcV) } ,
        { minDy(x, y, _tcU), minDy(x, y, _tcV) }
      };

      // ONLY compute the eigenvalues after checking that the matrix
      // is nonsingular
      JAMA::LU<float> LU = computeLU2x2(jacobian);

      if (LU.isNonsingular())
      {
        // get the analytic eigenvalues, quite slow right now...
        float eigenvalues[2] = {1., 1.};
        computeEigenvalues2x2(&eigenvalues[0], jacobian);
        _eigMax(x, y) = (eigenvalues[0] > eigenvalues[1]) ? eigenvalues[0] : eigenvalues[1];
        _eigMin(x, y) = (eigenvalues[0] < eigenvalues[1]) ? eigenvalues[0] : eigenvalues[1];
        maxeig = (maxeig > _eigMax(x, y)) ? maxeig : _eigMax(x, y);
        maxeig = (mineig < _eigMax(x, y)) ? mineig : _eigMax(x, y);
      }
      else
      {
        _eigMax[index] = 10.0f;
        _eigMin[index] = 0.1;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// advect & reset texture coordinates based on eigenvalues
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::resetTextureCoordinates()
{
  // allowed deformation of the textures
  const float limit = 2.f;
  const float limitInv = 1./limit;

  // standard reset
  int resets = 0;
  const float dx = 1./(float)(_resSm[0]);
  const float dy = 1./(float)(_resSm[1]);

  for (int y = 1; y < _yResSm - 1; y++)
    for (int x = 1; x < _xResSm - 1; x++)
    {
      if (_eigMax(x, y) > limit || _eigMin(x, y) < limitInv)
      {
        _tcU(x, y) = (float) x * dx;
        _tcV(x, y) = (float) y * dy;
        resets++;
      }
    }
}

//////////////////////////////////////////////////////////////////////
// Compute the highest frequency component of the wavelet
// decomposition
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::decomposeEnergy()
{
  // do the decomposition -- the goal here is to have
  // the energy with the high frequency component stomped out
  // stored in _tcTemp when it is done. _highFreqEnergy is only used
  // as an additional temp array

  // downsample input
  downsampleXNeumann(_highFreqEnergy.data(), _energy.data(), _xResSm, _yResSm);
  downsampleYNeumann(_tcTemp.data(), _highFreqEnergy.data(), _xResSm, _yResSm);

  // upsample input
  upsampleYNeumann(_highFreqEnergy.data(), _tcTemp.data(), _xResSm, _yResSm);
  upsampleXNeumann(_tcTemp.data(), _highFreqEnergy.data(), _xResSm, _yResSm);

  // subtract the down and upsampled field from the original field --
  // what should be left over is solely the high frequency component
  for (int y = 0; y < _yResSm; y++)
    for (int x = 0; x < _xResSm; x++) {
      // brute force reset of boundaries
      if(x >= _xResSm - 1 || y >= _yResSm - 1 || y <= 0 || x <= 0)
        _highFreqEnergy(x, y) = 0.;
      else
        _highFreqEnergy(x, y) = _energy(x, y) - _tcTemp(x, y);
  }
}

//////////////////////////////////////////////////////////////////////
// compute velocity from energies
// for wavelet decomposition
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::computeEnergy(FLUID_2D& xvel, FLUID_2D& yvel)
{
  // compute everywhere
  for (int y = 0; x < _yResSm; y++)
    for (int x = 0; x < _xResSm; x++)
      _energy[x] = 0.5f * (xvel(x, y) * xvel(x, y) + yvel(x, y) * yvel(x, y));

  FLUID_2D_BOUNDED::copyBoundary(_energy, true, true, 1);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Evaluate derivatives
//////////////////////////////////////////////////////////////////////////////////////////
void WTURBULENCE::WVelocity(float* orgPos, float* out)
{
  out[0] = WNoiseDy(orgPos, _noiseTile.data());
  out[1] = -WNoiseDx(orgPos, _noiseTile.data());
}

//////////////////////////////////////////////////////////////////////
// perform an actual noise advection step
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::stepTurbulenceReadable(float dtOrg, FLUID_2D& xvel, FLUID_2D& yvel)
{
  // enlarge timestep to match grid
  const float dt = dtOrg * _amplify;
  const float invAmp = 1.0f / _amplify;

  // prepare textures
  advectTextureCoordinates(dtOrg, xvel, yvel);

  // compute eigenvalues of the texture coordinates
  computeEigenvalues();

  // do wavelet decomposition of energy
  computeEnergy(xvel, yvel);
  decomposeEnergy();

  float maxVelocity = 0.;
  for (int y = 1; y < _yResBig - 1; y++)
    for (int x = 1; x < _xResBig - 1; x++)
    {
      // get unit position for both fine and coarse grid
      const int pos[2] = {x, y};
      const int posSm[2] = {(int) x * invAmp, (int) y * invAmp};

      // get a linearly interpolated velocity and texcoords
      // from the coarse grid
      float vel[2] = {};
      float uv[2] = {};
      INTERPOLATE::lerp2dVec(xvel, yvel, posSm[0], posSm[1], vel);
      INTERPOLATE::lerp3dVec(_tcU, _tcV, posSm[0], posSm[1], uv);

      // multiply the texture coordinate by _resSm so that turbulence
      // synthesis begins at the first octave that the coarse grid
      // cannot capture
      float texCoord[2] = {uv[0] * _resSm[0],  uv[1] * _resSm[1]};

      // retrieve wavelet energy at highest frequency
      float energy = INTERPOLATE::lerp2d(
          _highFreqEnergy, posSm[0], posSm[1]);

      // base amplitude for octave 0
      float coefficient = sqrtf(2.0f * fabs(energy));
      const float amplitude = _strength * fabs(0.5 * coefficient) * persistence;

      // add noise to velocity, but only if the turbulence is
      // sufficiently undeformed, and the energy is large enough
      // to make a difference
      const bool addNoise = _eigMax[indexSmall] < 2. &&
                            _eigMin[indexSmall] > 0.5;
      if (addNoise && amplitude > _cullingThreshold) {
        // base amplitude for octave 0
        float amplitudeScaled = amplitude;

        for (int octave = 0; octave < _octaves; octave++)
        {
          // multiply the vector noise times the maximum allowed
          // noise amplitude at this octave, and add it to the total
          float wvel[2] = {};
          WVelocity(texCoord, wvel);
          vel[0] += wvel[0] * amplitudeScaled;
          vel[1] += wvel[1] * amplitudeScaled;

          // scale coefficient for next octave
          amplitudeScaled *= persistence;
          texCoord[0] *= 2.0f;
          texCoord[1] *= 2.0f;
        }
      }

      // Store velocity + turbulence in big grid for maccormack step
      //
      // If you wanted to save memory, you would instead perform a
      // semi-Lagrangian backtrace for the current grid cell here. Then
      // you could just throw the velocity away.
      _bigUx(x, y) = vel[0];
      _bigUy(x, y) = vel[1];

      // compute the velocity magnitude for substepping later
      const float velMag = _bigUx(x, y) * _bigUx(x, y) + _bigUy(x, y) * _bigUy(x, y);
      if (velMag > maxVelocity) maxVelocity = velMag;
    }

  // prepare density for an advection
  FLUID_2D_BOUNDED::swapFields(_densityBig, _densityBigOld);

  // based on the maximum velocity present, see if we need to substep,
  // but cap the maximum number of substeps to 5
  const int maxSubSteps = 5;
  maxVelocity = sqrt(maxVelocity) * dt;
  int totalSubsteps = (int)(maxVelocity / (float)maxSubSteps);
  totalSubsteps = (totalSubsteps < 1) ? 1 : totalSubsteps;
  totalSubsteps = (totalSubsteps > maxSubSteps) ? maxSubSteps : totalSubsteps;
  const float dtSubdiv = dt / (float)totalSubsteps;

  // set boundaries of big velocity grid
  FLUID_2D_BOUNDED::fillBoundary(_bigUx, true, false);
  FLUID_2D_BOUNDED::fillBoundary(_bigUy, false, true);

  // do advection, with substepping if necessary
  for(int substep = 0; substep < totalSubsteps; substep++)
  {
    FLUID_2D_BOUNDED::advect(dtSubdiv, _densityBig, _densityBigOld, _bigUx, _bigUy);

    if (substep < totalSubsteps - 1)
      FLUID_2D_BOUNDED::swapFields(_densityBig, _densityBigOld);
  } // substep

  // wipe the density borders
  FLUID_3D::fillBoundary(_densityBig);

  // reset texture coordinates now in preparation for next timestep
  // Shouldn't do this before generating the noise because then the
  // eigenvalues stored do not reflect the underlying texture coordinates
  resetTextureCoordinates();

  _totalStepsBig++;
}

