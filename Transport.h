#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "TRandom3.h"
#include "OpticsConfig.h"
#include "Math/Vector3D.h"
#include <vector>

using V3 = ROOT::Math::XYZVector;

PhotonResult PropagateOnePhoton(
    TRandom3 &rng,
    const V3 &sitePos,
    int site_number,
    const OpticsConfig &cfg,
    const std::vector<V3> &normals,
    const std::vector<V3> &pointPlane
);

// Finds the closest valid plane hit.
// Returns plane index, or -1 if none, or 999 if multiple (your current logic).
int getIntersectionPlan(
    const std::vector<V3> &normals,
    const std::vector<V3> &pointPlane,
    const V3 &start_point,
    const V3 &direction,
    double &tOut,
    V3 &pOut,
    const OpticsConfig &cfg,
    double eps = 1e-12
);

// Plane intersection primitive
bool getIntersectionPoint(
    const V3 &planeNormal,
    const V3 &planePoint,
    const V3 &start_point,
    const V3 &direction,
    double &tOut,
    V3 &pOut,
    double eps = 1e-12
);

#endif
