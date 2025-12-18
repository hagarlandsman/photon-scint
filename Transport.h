#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "TRandom3.h"
#include "Vec3.h"
#include "OpticsConfig.h"

PhotonResult PropagateOnePhoton(TRandom3 &rng, const Vec3 &sitePos, int site_number, const OpticsConfig &cfg);

#endif
