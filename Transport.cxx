#include "Transport.h"
#include "TRandom3.h"
#include "Vec3.h"
#include "OpticsConfig.h"
#include "Geometry.h"
#include "TMath.h"
#include <cmath>
#include <algorithm>

PhotonResult PropagateOnePhoton(TRandom3 &rng, const Vec3 &sitePos, int site_number, const OpticsConfig &cfg);

PhotonResult PropagateOnePhoton(
    TRandom3 &rng,
    const Vec3 &sitePos,
    int site_number,
    const OpticsConfig &cfg)
{
    PhotonResult r;
    r.site_number = site_number;
    r.x0 = sitePos.x;
    r.y0 = sitePos.y;
    r.z0 = sitePos.z;
    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;

    const int PTFE = (wrap == 1);
    const int mylar = (wrap == 2);

    double sinThetaC = cfg.nOut / cfg.nScint;
    if (sinThetaC >= 1.0)
        sinThetaC = 0.999999;
    double cosThetaC = std::cos(std::asin(sinThetaC));

    Vec3 pos = sitePos;
    Vec3 dir = Unit(SampleIsotropic(rng));

    if (cfg.savePath)
    {
        r.xPath.clear();
        r.yPath.clear();
        r.zPath.clear();
        r.xPath.push_back((float)pos.x);
        r.yPath.push_back((float)pos.y);
        r.zPath.push_back((float)pos.z);
    }

    for (int step = 0; step < cfg.maxSteps; step++)
    {

        // Candidate intersections with standard box
        double tx = 1e99;
        int xPlane = -1;
        if (dir.x > 1e-15)
        {
            tx = (L - pos.x) / dir.x;
            xPlane = 1;
        }
        if (dir.x < -1e-15)
        {
            tx = (0.0 - pos.x) / dir.x;
            xPlane = 0;
        }

        double tyP = 1e99, tyM = 1e99;
        if (dir.y > 1e-15)
            tyP = (W / 2.0 - pos.y) / dir.y;
        if (dir.y < -1e-15)
            tyM = (-W / 2.0 - pos.y) / dir.y;

        double tzP = 1e99, tzM = 1e99;
        if (dir.z > 1e-15)
            tzP = (T / 2.0 - pos.z) / dir.z;
        if (dir.z < -1e-15)
            tzM = (-T / 2.0 - pos.z) / dir.z;

        // Optional wedge planes: y = +/- yMax(x) in wedge regions (taper in y only)
        double tWedgeP = 1e99, tWedgeM = 1e99;
        int wedgeSide = 0; // 0 none, 1 left, 2 right

        if (cfg.useWedge && cfg.wedgeLen > 0 && cfg.wedgeTipW > 0 && cfg.wedgeTipW <= W)
        {

            if (InLeftWedge(pos.x, cfg.wedgeLen))
                wedgeSide = 1;
            else if (InRightWedge(pos.x, L, cfg.wedgeLen))
                wedgeSide = 2;

            if (wedgeSide != 0)
            {

                double a = 0.5 * cfg.wedgeTipW;
                double b = (0.5 * W - 0.5 * cfg.wedgeTipW) / cfg.wedgeLen; // >= 0

                // For left: yMax = a + b*x
                // For right: yMax = a + b*(L-x) = (a + b*L) - b*x
                // We solve intersections with:
                //   y(t) = +yMax(x(t))
                //   y(t) = -yMax(x(t))

                // + surface: y - yMax(x) = 0
                {
                    double denom = 0.0;
                    double rhs = 0.0;

                    if (wedgeSide == 1)
                    {
                        // y - (a + b x) = 0 -> t*(dy - b*dx) = a + b*x0 - y0
                        denom = dir.y - b * dir.x;
                        rhs = a + b * pos.x - pos.y;
                    }
                    else
                    {
                        // y - ((a + bL) - b x) = 0 -> y - a0 + b x = 0
                        // t*(dy + b*dx) = a0 - b*x0 - y0, with a0 = a + bL
                        double a0 = a + b * L;
                        denom = dir.y + b * dir.x;
                        rhs = a0 - b * pos.x - pos.y;
                    }

                    if (std::fabs(denom) > 1e-15)
                    {
                        double t = rhs / denom;
                        if (t > 1e-12)
                        {
                            double xh = pos.x + t * dir.x;
                            double zh = pos.z + t * dir.z;
                            bool okx = (wedgeSide == 1) ? (xh >= 0.0 && xh <= cfg.wedgeLen)
                                                        : (xh >= L - cfg.wedgeLen && xh <= L);
                            if (okx && std::fabs(zh) <= T / 2.0)
                                tWedgeP = t;
                        }
                    }
                }

                // - surface: y + yMax(x) = 0
                {
                    double denom = 0.0;
                    double rhs = 0.0;

                    if (wedgeSide == 1)
                    {
                        // y + (a + b x) = 0 -> t*(dy + b*dx) = -(a + b*x0) - y0
                        denom = dir.y + b * dir.x;
                        rhs = -(a + b * pos.x) - pos.y;
                    }
                    else
                    {
                        // y + ((a + bL) - b x) = 0 -> y + a0 - b x = 0
                        // t*(dy - b*dx) = -a0 + b*x0 - y0
                        double a0 = a + b * L;
                        denom = dir.y - b * dir.x;
                        rhs = -a0 + b * pos.x - pos.y;
                    }

                    if (std::fabs(denom) > 1e-15)
                    {
                        double t = rhs / denom;
                        if (t > 1e-12)
                        {
                            double xh = pos.x + t * dir.x;
                            double zh = pos.z + t * dir.z;
                            bool okx = (wedgeSide == 1) ? (xh >= 0.0 && xh <= cfg.wedgeLen)
                                                        : (xh >= L - cfg.wedgeLen && xh <= L);
                            if (okx && std::fabs(zh) <= T / 2.0)
                                tWedgeM = t;
                        }
                    }
                }

                // Inside wedge region, the wedge replaces the full y=+-W/2 planes
                tyP = 1e99;
                tyM = 1e99;
            }
        }

        // Pick the nearest intersection
        double tmin = 1e99;
        int plane = -1;
        if (tx > 1e-12 && tx < tmin)
        {
            tmin = tx;
            plane = xPlane;
        }
        if (tyM > 1e-12 && tyM < tmin)
        {
            tmin = tyM;
            plane = 2;
        }
        if (tyP > 1e-12 && tyP < tmin)
        {
            tmin = tyP;
            plane = 3;
        }
        if (tzM > 1e-12 && tzM < tmin)
        {
            tmin = tzM;
            plane = 4;
        }
        if (tzP > 1e-12 && tzP < tmin)
        {
            tmin = tzP;
            plane = 5;
        }
        if (tWedgeP > 1e-12 && tWedgeP < tmin)
        {
            tmin = tWedgeP;
            plane = 6;
        }
        if (tWedgeM > 1e-12 && tWedgeM < tmin)
        {
            tmin = tWedgeM;
            plane = 7;
        }

        if (plane < 0 || tmin > 1e98)
            break;

        // Bulk absorption
        if (cfg.absLen > 0)
        {
            double survive = std::exp(-tmin / cfg.absLen);
            if (rng.Uniform() > survive)
            {
                r.path += tmin;
                r.absorbed = 1;
                r.endPlane = plane;

                r.xf = pos.x + dir.x * tmin;
                r.yf = pos.y + dir.y * tmin;
                r.zf = pos.z + dir.z * tmin;

                if (cfg.savePath)
                {
                    r.xPath.push_back((float)r.xf);
                    r.yPath.push_back((float)r.yf);
                    r.zPath.push_back((float)r.zf);
                }
                return r;
            }
        }

        // Move to boundary
        pos = pos + dir * tmin;
        r.path += tmin;
        r.endPlane = plane;
        r.xf = pos.x;
        r.yf = pos.y;
        r.zf = pos.z;
        if (cfg.savePath)
        {
            r.xPath.push_back((float)pos.x);
            r.yPath.push_back((float)pos.y);
            r.zPath.push_back((float)pos.z);
        }
        if (!InsideActiveVolume(pos.x, pos.y, pos.z, cfg))
        {
            r.escaped = 1;
            r.endPlane = plane;
            return r;
        }

        // End planes: detection
        if (plane == 0 || plane == 1)
        {
            if (!InsideWedgeAperture(pos.x, pos.y, pos.z, cfg))
            {
                r.escaped = 1;
                r.endPlane = plane;
                return r;
            }

            r.reachedEnd = 1;
            r.pmt_side = (plane == 0) ? 1 : 2;
            r.inPMT = InsidePMTCircle(pos.y, pos.z, cfg.rPMT) ? 1 : 0;

            r.epsRaysPMT = epsCoupleExp(pos.y, pos.z, cfg.rPMT, cfg.eps0, cfg.lambdaC);
            double pEnd = cfg.epsCouple * cfg.pde * r.epsRaysPMT;
            if (rng.Uniform() < pEnd)
                r.detected = 1;
          //  printf("Detected\n");
            return r;
        }

        // Surface normal
        Vec3 nHat;
        if (plane == 2)
            nHat = Vec3(0, -1, 0);
        if (plane == 3)
            nHat = Vec3(0, +1, 0);
        if (plane == 4)
            nHat = Vec3(0, 0, -1);
        if (plane == 5)
            nHat = Vec3(0, 0, +1);

        if (plane == 6 || plane == 7)
        {
            double b = (cfg.wedgeLen > 0) ? (0.5 * W - 0.5 * cfg.wedgeTipW) / cfg.wedgeLen : 0.0;

            if (plane == 6)
            {
                // + surface: y = +yMax(x)
                // Left: F = y - a - b x -> grad = (-b, +1, 0)
                // Right: F = y - (a + bL) + b x -> grad = (+b, +1, 0)
                double nx = (wedgeSide == 2) ? +b : -b;
                nHat = Unit(Vec3(nx, +1.0, 0.0));
            }
            else
            {
                // - surface: y = -yMax(x)
                // Left: F = y + a + b x -> outward should be toward -y
                // Use normal pointing outward: (+b, -1, 0)
                // Right: F = y + (a + bL) - b x -> outward: (-b, -1, 0)
                double nx = (wedgeSide == 2) ? -b : +b;
                nHat = Unit(Vec3(nx, -1.0, 0.0));
            }
        }

        nHat = Unit(nHat);

        double cosInc = std::fabs(Dot(dir, nHat));
        bool isTIR = (cosInc < cosThetaC);

        if (isTIR)
        {
            dir = Unit(ReflectSpecular(dir, nHat));
            r.nBounces++;
            continue;
        }

        // Non-TIR: reflect with probability Rwrap, else escape
        if (rng.Uniform() >= cfg.Rwrap)
        {
            r.escaped = 1;
            return r;
        }

        // Wrap model: PTFE -> diffuse, Mylar -> specular
        if (PTFE && !mylar)
            dir = SampleLambert(-nHat, rng);
        else
            dir = Unit(ReflectSpecular(dir, nHat));

        r.nBounces++;
    }

    return r;
}
