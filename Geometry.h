#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "OpticsConfig.h"
bool InsideActiveVolume(double x, double y, double z, const OpticsConfig &cfg);
bool InsideWedgeAperture(double x, double y, double z, const OpticsConfig &cfg);

bool InsidePMTCircle(double y, double z, double rPMT);
double epsCoupleExp(double y, double z, double rPMT, double eps0, double lambdaC);

bool InLeftWedge(double x, double Lg);
bool InRightWedge(double x, double L, double Lg);



#endif
