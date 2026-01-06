#include "Geometry.h"
#include <cmath>
#include <algorithm>
#include "Math/Vector3D.h"

using V3 = ROOT::Math::XYZVector;

 bool InsideActiveVolume(double x, double y, double z,
                                      const OpticsConfig &cfg)
{

    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;
    double eps = 1e-9;
    if (z < -T / 2.0  - eps|| z > T / 2.0 + eps ){
        printf (" z=%f out of bounds (%f, %f)\n",z,-T/2.0,T/2.0);
        return false;
    }

    if (!cfg.useWedge || cfg.wedgeLen <= 0)
    {
        return (y >= -W / 2.0 - eps && y <= W / 2.0+eps) && (x >= 0.0 -eps && x <= L+eps);
    }

    if (x > 0.5 * L + eps || x < -0.5*L-eps) {
        printf (" x=%f out of bounds (%f, %f)\n",x,0.0,L);
        return false;
    }

    // left wedge [0, wedgeLen]
    if (x <= cfg.wedgeLen + eps)
    {
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * (x / cfg.wedgeLen);
        bool ok = std::fabs(y) <= yMax;
        if (not ok)
            printf (" y=%f out of bounds (%f, %f) in left wedge\n",y,-yMax,yMax);
        return ok ;
    }

    // right wedge [L-wedgeLen, L]
    if (x >= L - cfg.wedgeLen-eps  )
    {
        double u = (L - x) / cfg.wedgeLen;
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * u;
       bool ok = std::fabs(y) <= yMax;
         if (not ok)
            printf (" y=%f out of bounds (%f, %f) in right wedge\n",y,-yMax,yMax);
        return ok;
    }

    // middle region
    return std::fabs(y) <= W / 2.0;
}

bool InsidePMTCircle(double y, double z, double rPMT)
{
    return (y * y + z * z) <= (rPMT * rPMT);
}

 double epsCoupleExp(double y, double z, double rPMT, double eps0, double lambdaC)
{
    double d = std::sqrt(y * y + z * z);
    if (d <= rPMT)
        return 1.0;
    return eps0 * std::exp(-d / lambdaC);
}

// Wedge logic: taper in y only near both ends.
// Left wedge active for x in [0, Lg], right wedge for x in [L-Lg, L].
// Width shrinks linearly from W at x=Lg (or x=L-Lg) to wTip at x=0 (or x=L).

 bool InLeftWedge(double x, double Lg) { return x >= 0.0 && x <= Lg; }
 bool InRightWedge(double x, double L, double Lg) { return x >= (L - Lg) && x <= L; }

 double WedgeYMax(double x, double L, double W, double Lg, double wTip, int side)
{
    // side: 1 = left wedge, 2 = right wedge
    double a = 0.5 * wTip;
    double b = (0.5 * W - 0.5 * wTip) / Lg; // >= 0

    if (side == 1)
    {
        // x in [0, Lg]
        return a + b * x;
    }
    else
    {
        // x in [L-Lg, L]
        return a + b * (L - x);
    }
}
 bool InsideWedgeAperture(double x, double y, double z,
                                       const OpticsConfig &cfg)
{
    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;
    const double eps = 1e-6;
    if (std::fabs(z) > T / 2.0 + eps)
        return false;
    if (x < -L *0.5 - eps || x > L*0.5 + eps)
        return false;

    if (!cfg.useWedge || cfg.wedgeLen <= 0)
    {
        return std::fabs(y) <= W / 2.0 + eps;
    }
    if (cfg.wedgeTipW <= 0.0 || cfg.wedgeTipW > W)
    {
        // treat as no wedge (or return false, depending on your preference)
        return std::fabs(y) <= W / 2.0 + eps;
    }

    // Left wedge: x in [0, wedgeLen]
    if (x <= cfg.wedgeLen + eps)
    {
        double u = std::max(0.0, std::min(1.0, x / cfg.wedgeLen));
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * u;
        return std::fabs(y) <= yMax + eps;
    }
    // Right wedge: x in [L-wedgeLen, L]
    if (x >= L - cfg.wedgeLen - eps)
    {
        double u = std::max(0.0, std::min(1.0, (L - x) / cfg.wedgeLen));
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * u;
        return std::fabs(y) <= yMax + eps;
    }

    // middle rectangular region
    return std::fabs(y) <= W / 2.0 + eps;
}

 V3 ReflectSpecular(const V3 v, const V3 nHat)
{
    return v - nHat * (2.0 * v.Dot(nHat));
}

 V3 SampleIsotropic(TRandom3 &rng)
{
    double u = rng.Uniform(-1.0, 1.0);
    double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    double s = std::sqrt(std::max(0.0, 1.0 - u * u));
    return V3(s * std::cos(phi), s * std::sin(phi), u);
}

// Cosine-weighted (Lambertian) hemisphere around nHat

 V3 SampleLambert(const V3 &nHat, TRandom3 &rng)
{
    V3 w = nHat.Unit();
    V3 a = (std::fabs(w.x()) < 0.9) ? V3(1, 0, 0) : V3(0, 1, 0);
    V3 u = V3(
        w.y() * a.z() - w.z() * a.y(),
        w.z() * a.x() - w.x() * a.z(),
        w.x() * a.y() - w.y() * a.x()).Unit();
    V3 v = V3(
        w.y() * u.z() - w.z() * u.y(),
        w.z() * u.x() - w.x() * u.z(),
        w.x() * u.y() - w.y() * u.x());

    double r1 = rng.Uniform();
    double r2 = rng.Uniform();
    double phi = 2.0 * TMath::Pi() * r1;
    double cosTheta = std::sqrt(1.0 - r2);
    double sinTheta = std::sqrt(r2);

    V3 dirLocal(
        std::cos(phi) * sinTheta,
        std::sin(phi) * sinTheta,
        cosTheta);

    V3 dir = u * dirLocal.x() + v * dirLocal.y() + w * dirLocal.z();
    return dir.Unit();
}