#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "OpticsConfig.h"
#include "Math/Vector3D.h"
#include "TRandom3.h"
#include "TMath.h"

using V3 = ROOT::Math::XYZVector;

bool InsideActiveVolume(double x, double y, double z, const OpticsConfig &cfg);
bool InsideWedgeAperture(double x, double y, double z, const OpticsConfig &cfg);

bool InsidePMTCircle(double y, double z, double rPMT);
double epsCoupleExp(double y, double z, double rPMT, double eps0, double lambdaC);

bool InLeftWedge(double x, double Lg);
bool InRightWedge(double x, double L, double Lg);

 V3 ReflectSpecular(const V3 v, const V3 nHat);

 V3 SampleIsotropic(TRandom3 &rng);
V3 SampleLambert(const V3 &nHat, TRandom3 &rng);

#endif
