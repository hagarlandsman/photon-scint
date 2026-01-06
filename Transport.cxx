#include "Transport.h"
#include "Geometry.h"
#include "TMath.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using V3 = ROOT::Math::XYZVector;

// You call this, so define it here (or include the header that defines it).


PhotonResult PropagateOnePhoton(
    TRandom3 &rng,
    const V3 &sitePos,
    int site_number,
    const OpticsConfig &cfg,
    const std::vector<V3> &normals,
    const std::vector<V3> &pointPlane
)
{

    bool debug=false;
    if (debug) printf ("  sitePos=(%f, %f, %f), site_number=%d\n", sitePos.x(), sitePos.y(), sitePos.z(), site_number);
    PhotonResult r;
    r.site_number = site_number;
    r.x0 = sitePos.x();
    r.y0 = sitePos.y();
    r.z0 = sitePos.z();

    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;

    double sinThetaC = cfg.nOut / cfg.nScint;
    if (sinThetaC >= 1.0) sinThetaC = 0.999999;
    double cosThetaC = std::cos(std::asin(sinThetaC));

    V3 pos = sitePos;
    V3 dir = SampleIsotropic(rng).Unit();
    //dir=V3(0,1,0);
    if (cfg.savePath)
    {
        if (debug)  cout<<"Saving path...\n";
        r.xPath.clear();
        r.yPath.clear();
        r.zPath.clear();
        r.xPath.push_back((float)pos.x());
        r.yPath.push_back((float)pos.y());
        r.zPath.push_back((float)pos.z());
    }
    for (int step = 0; step < cfg.maxSteps; ++step)
    {
        double tmin = 0.0;
        V3 hitPoint;
       if (debug)  cout<<"Step "<<step<<", pos=("<<pos.x()<<", "<<pos.y()<<", "<<pos.z()<<"), dir=("<<dir.x()<<", "<<dir.y()<<", "<<dir.z()<<")\n";
        int plane = getIntersectionPlan(normals, pointPlane, pos, dir, tmin, hitPoint, cfg);
        if (plane >= 0)
          if (debug)   cout<<"\tHit plane "<<plane<<" at tmin="<<tmin<<", point=("<<hitPoint.x()<<", "<<hitPoint.y()<<", "<<hitPoint.z()<<")\n";
        if (plane < 0) {
           if (debug)  cout<<"\t No intersection found, photon escapes. XXXX\n";
            break;
        }
        pos = hitPoint;

        if (cfg.savePath)
                {
                    r.xPath.push_back((float)hitPoint.x());
                    r.yPath.push_back((float)hitPoint.y());
                    r.zPath.push_back((float)hitPoint.z());
                }
        // Bulk absorption
        if (cfg.absLen > 0)
        {
            double survive = std::exp(-tmin / cfg.absLen);
            if (rng.Uniform() > survive)
            {
                r.path += tmin;
                r.absorbed = 1;
                r.endPlane = plane;

                r.xf = hitPoint.x();
                r.yf = hitPoint.y();
                r.zf = hitPoint.z();


               if (debug)  cout<<"\t\tPhoton absorbed in bulk at step "<<step<<", plane "<<plane<<"\n";
                return r;
            }
        }

        // Move
//        pos = pos + tmin * dir;
        r.path += tmin;
        r.endPlane = plane;

        r.xf = pos.x();
        r.yf = pos.y();
        r.zf = pos.z();

        if (!isInside(pos, cfg))
        {
            r.escaped = 1;
            if (debug)  cout<<"\t\tPhoton escaped active volume at step "<<step<<", plane "<<plane<<"\n";
            return r;
        }

        // End planes: detection
        if (plane == 0 || plane == 1)
        {
            /* if (!InsideWedgeAperture(pos.x(), pos.y(), pos.z(), cfg))

            {
                r.escaped = 1;
                cout<<"\t\tPhoton escaped through wedge aperture at step "<<step<<", plane "<<plane<<"\n";
                return r;
            }
            */
            r.inPMT = InsidePMTCircle(pos.y(), pos.z(), cfg.rPMT) ? 1 : 0;
            if (debug) printf ("in plane0,1 z=%f \n",pos.z());
            if (r.inPMT)
            {
                if (debug) cout<<"\t\tPhoton inside PMT circle at step "<<step<<", plane "<<plane<<"\n";
                r.reachedEnd = 1;
                r.pmt_side = (plane == 0) ? 1 : 2;

                r.epsRaysPMT = epsCoupleExp(pos.y(), pos.z(), cfg.rPMT, cfg.eps0, cfg.lambdaC);
                double pEnd = cfg.epsCouple * cfg.pde * r.epsRaysPMT;

                r.detected = (rng.Uniform() < pEnd) ? 1 : 0;
                if (debug)  cout<<"\t\tPhoton detected at step "<<step<<", plane "<<plane<<"\n";
                return r;
            }
        }

        // Reflection
        V3 nHat = normals[plane].Unit();
        double cosInc = std::fabs(dir.Dot(nHat));
        bool isTIR = (cosInc < cosThetaC);

        if (isTIR)
        {
           if (debug)  cout<<"\t\tTotal internal reflection at step "<<step<<", plane "<<plane<<"\n";
            dir = ReflectSpecular(dir, nHat).Unit();
            r.nBounces++;
            continue;
        }
       if (debug)  cout<<"\t\tNot TIR at step "<<step<<", plane "<<plane<<"\n";
        if (rng.Uniform() >= cfg.Rwrap )
        {
          if (debug)   cout<<"\t\tPhoton escaped lost reflection at step "<<step<<", plane "<<plane<<"  Rwrap="<<cfg.Rwrap<<"\n";
            r.escaped = 1;
            return r;
        }

        if (cfg.wrap == 1){
            dir = SampleLambert(nHat, rng).Unit();
          if (debug)   cout<<"\t\tDiffuse reflection (PTFE) at step "<<step<<", plane "<<plane<<"\n";
                }
        else
            {dir = ReflectSpecular(dir, nHat).Unit();
           if (debug)  cout<<"Specular reflection (Mylar) at step "<<step<<", plane "<<plane<<"\n";}
        r.nBounces++;
    }
    cout<<"Warning: max steps reached in photon propagation"<<cfg.maxSteps<<"\n";
    return r;
}

/*
int getIntersectionPlan(
    const std::vector<V3> &normals,
    const std::vector<V3> &pointPlane,
    const V3 &start_point,
    const V3 &direction,
    double &tOut,
    V3 &pOut,
    const OpticsConfig &cfg,
    double eps
)
{
    int bestPlane = -1;
    double bestT = 1e99;
    V3 bestP(0,0,0);

    int nValid = 0;

    int N = (int)normals.size();
    for (int i = 0; i < N; ++i)
    {
        double t;
        V3 hitPoint;

        bool ok = getIntersectionPoint(normals[i], pointPlane[i], start_point, direction, t, hitPoint, eps);
        if (!ok) continue;

        // Must be inside the finite face polygon/extent (your isInside)
        int loc = isInside(hitPoint, cfg.L, cfg.W, cfg.T, cfg.useWedge ? cfg.wedgeLen : 0.0, cfg.useWedge ? cfg.wedgeTipW : 0.0);

        if (loc > 0)
        {gr
            nValid++;

            if (t < bestT)
            {
                bestT = t;
                bestPlane = i;
                bestP = hitPoint;
            }
        }
    }

    if (bestPlane < 0)
        return -1;

    tOut = bestT;
    pOut = bestP;

    // if you still want to flag "multi hit" cases:
    if (nValid > 1)
        return 999;

    return bestPlane;
}

*/

inline int getIntersectionPlan(const std::vector<V3> &normals,
                              const std::vector<V3> &pointPlane,
                              const V3 &start_point,
                              const V3 &direction,
                              double &tOut,
                              V3 &pOut,
                              const OpticsConfig &cfg,
                              double eps )
{
    int bestPlane = -1;
    double bestT = 1e99;
    V3 bestP;

    for (int i = 0; i < (int)normals.size(); ++i)
    {
        double t = 0.0;
        V3 hitPoint;
        bool ok = getIntersectionPoint(normals[i], pointPlane[i],
                                       start_point, direction,
                                       t, hitPoint, eps);

      //  cout<<" Checking plane "<<i<<": t="<<t<<"  hitPoint=("<<hitPoint.x()<<","<<hitPoint.y()<<","<<hitPoint.z()<<") \n";
        if (!ok)
            continue;

        if (t <= eps)
            continue;


        //int loc = isInside(hitPoint, cfg.L, cfg.W, cfg.T, cfg.wedgeLen, cfg.wedgeTipW, 1e-12);
       // if (loc <= 0) continue;

        if (t < bestT)
        {
            bestT = t;
            bestPlane = i;
            bestP = hitPoint;
        }
    }

    if (bestPlane < 0) return -1;

    tOut = bestT;
    pOut = bestP;
    return bestPlane;
}


bool getIntersectionPoint(
    const V3 &planeNormal,
    const V3 &planePoint,
    const V3 &start_point,
    const V3 &direction,
    double &tOut,
    V3 &pOut,
    double eps
)
{
    tOut = 0.0;

    double denom = planeNormal.Dot(direction);
    if (std::fabs(denom) < eps)
        return false;

    double t = planeNormal.Dot(planePoint - start_point) / denom;
    if (t < 0)
        return false;

    tOut = t;
    pOut = start_point + t * direction;
    return true;
}
