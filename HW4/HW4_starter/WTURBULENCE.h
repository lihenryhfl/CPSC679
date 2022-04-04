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

#ifndef WTURBULENCE_H
#define WTURBULENCE_H

#include "FIELD_2D.h"
using namespace BasicVector;
class SIMPLE_PARSER;

///////////////////////////////////////////////////////////////////////////////
/// Main WTURBULENCE class, stores large density array etc.
///////////////////////////////////////////////////////////////////////////////
class WTURBULENCE
{
  public:
    // both config files can be NULL, altCfg might override values from noiseCfg
    WTURBULENCE(int xResSm, int yResSm, int zResSm, int amplify);

    /// destructor
    virtual ~WTURBULENCE();

    // step more readable version -- no rotation correction
    void stepTurbulenceReadable(float dt, float* xvel, float* yvel, float* zvel, unsigned char *obstacles);

    // texcoord functions
    void advectTextureCoordinates(float dtOrg, float* xvel, float* yvel, float* zvel);
    void resetTextureCoordinates();

    void computeEnergy(float* xvel, float* yvel, float* zvel, unsigned char *obstacles);

    // evaluate wavelet noise function
    void WVelocity(float* orgPos, float* out);

    // access functions
    inline FIELD_2D* getDensityBig() { return &_densityBig; }
    inline FIELD_2D* getArrayTcU() { return &_tcU; }
    inline FIELD_2D* getArrayTcV() { return &_tcV; }
    inline FIELD_2D* getArrayTcW() { return &_tcW; }
    inline FIELD_2D* getArrayEigMin() { return &_eigMin; }
    inline FIELD_2D* getArrayEigMax() { return &_eigMax; }

    inline FIELD_2D* getResSm() { return &_resSm; }
    inline FIELD_2D* getResBig() { return &_resBig; }
    inline int getOctaves() { return _octaves; }

  protected:
    // enlargement factor from original velocity field / simulation
    // _Big = _amplify * _Sm
    int _amplify;
    int _octaves;
    float _strength;

    // noise settings
    float _cullingThreshold;
    float _noiseStrength;
    float _noiseSizeScale;
    bool _uvwAdvection;
    bool _uvwReset;
    float _noiseTimeanimSpeed;
    int _noiseControlType;
    // debug, scale density for projections output images
    float _outputScale;

    // noise resolution
    int _xResBig;
    int _yResBig;
    // original / small resolution
    int _xResSm;
    int _yResSm;
    float[3] _resSm;

    FIELD_2D _densityBig;
    FIELD_2D _densityBigOld;

    // big velocity macCormack fields
    FIELD_2D _bigUx;
    FIELD_2D _bigUy;
    // temp arrays for BFECC and MacCormack - they have more convenient
    // names in the actual implementations

    // texture coordinates for noise
    FIELD_2D _tcU;
    FIELD_2D _tcV;
    FIELD_2D _tcTemp;

    FIELD_2D _eigMin;
    FIELD_2D _eigMax;

    // wavelet decomposition of velocity energies
    FIELD_2D _energy;

    // noise data
    FIELD_2D _noiseTile;
    //float* _noiseTileExt;

    // step counter
    int _totalStepsBig;

    // highest frequency component of wavelet decomposition
    FIELD_2D _highFreqEnergy;

    void computeEigenvalues();
    void decomposeEnergy();
};

#endif // WTURBULENCE_H

