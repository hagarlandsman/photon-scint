#include "Geometry.h"
#include <cmath>
#include <algorithm>

 bool InsideActiveVolume(double x, double y, double z,
                                      const OpticsConfig &cfg)
{

    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;

    if (z < -T / 2.0 || z > T / 2.0)
        return false;

    if (!cfg.useWedge || cfg.wedgeLen <= 0)
    {
        return (y >= -W / 2.0 && y <= W / 2.0) && (x >= 0.0 && x <= L);
    }

    if (x < 0.0 || x > L)
        return false;

    // left wedge [0, wedgeLen]
    if (x <= cfg.wedgeLen)
    {
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * (x / cfg.wedgeLen);
        return std::fabs(y) <= yMax;
    }

    // right wedge [L-wedgeLen, L]
    if (x >= L - cfg.wedgeLen)
    {
        double u = (L - x) / cfg.wedgeLen;
        double yMax = 0.5 * cfg.wedgeTipW + (0.5 * W - 0.5 * cfg.wedgeTipW) * u;
        return std::fabs(y) <= yMax;
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
    const double eps = 1e-9;
    if (std::fabs(z) > T / 2.0 + eps)
        return false;
    if (x < 0.0 - eps || x > L + eps)
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

