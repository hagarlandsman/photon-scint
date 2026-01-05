#include <vector>
#include <iostream>
#include <cmath>

#include "TCanvas.h"
#include "TGraph.h"
#include "Math/Vector3D.h"

using V3 = ROOT::Math::XYZVector;

// Normal from 3 points: (p1-p0) x (p2-p0), normalized
static V3 findNormalFrom3Points(const V3 &p0, const V3 &p1, const V3 &p2)
{
    V3 n = (p1 - p0).Cross(p2 - p0);
    double mag = std::sqrt(n.Dot(n));
    if (mag > 0)
        n = (1.0 / mag) * n;
    double epsilon = 1e-15;
    if (std::fabs(n.x()) < epsilon)
        n.SetX(0.0);
    if (std::fabs(n.y()) < epsilon)
        n.SetY(0.0);
    if (std::fabs(n.z()) < epsilon)
        n.SetZ(0.0);
    // std::cout << "n = (" << n.x() << ", " << n.y() << ", " << n.z() << ")\n";
    return n;
}

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double eps = 1e-15);

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double wedgeLen, double wedgeTipW,
                    double eps = 1e-15);

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double eps)
{
    return isInside(p, L, W, T, 0.0, 0.0, eps);
}

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double wedgeLen, double wedgeTipW,
                    double eps)
{
    // z check
    if (p.z() < -0.5 * T - eps || p.z() > 0.5 * T + eps)
    {
        //  cout<<"z check failed\n";
        return 0;
    }
    // x check
    if (p.x() < 0.0 - eps || p.x() > L + 2 * wedgeLen + eps)
    {
        // cout<<"x check failed\n";
        return 0;
    }

    // y check
    if (p.y() < -0.5 * W - eps || p.y() > 0.5 * W + eps)
    {
        return 0;
    }

    // Main block return 1
    double x0 = wedgeLen;
    double x1 = L + wedgeLen;
    if (p.x() < x1 + eps && p.x() > x0 - eps && fabs(p.y()) <= 0.5 * W + eps && fabs(p.z()) <= 0.5 * T + eps)

        return 1;

    //    cout<<"Not in main block \n";
    double yMid = 0.5 * W;
    double yTip = 0.5 * wedgeTipW;

    // Left wedge (from (0,yTip) to (wedgeLen,yMid) )
    if (p.x() >= 0.0 - eps && p.x() <= wedgeLen + eps)
    {

        //    cout<<"in left wedge \n";

        double yMax = yTip + (yMid - yTip) * (p.x() / wedgeLen);

        if (std::fabs(p.y()) <= yMax + eps)
            return 2;
        else
            return 0;
    }

    // Right wedge (from (L+wedgeLen,yMid) to (L +2* wedgeLen,yTip) )
    if (p.x() >= L + wedgeLen - eps && p.x() <= L + 2 * wedgeLen + eps)
    {
        //  cout<<"in right wedge \n";
        double slope = (yMid - yTip) / (-wedgeLen);
        double yMax = yMid + slope * (p.x() - L - wedgeLen);
        cout<<"\n\nx= "<<p.x()<<"\tyMax= "<<yMax<<"\tslope= "<<slope<<"\tL0"<<L+wedgeLen-eps<<"\tL1="<<L+2*wedgeLen+eps<<"\n";
        if (std::fabs(p.y()) <= yMax + eps)
            return 2;
        else
            return 0;
    }

    //  cout<<"Not in wedge either \n";
    return 0;
}

#include <cmath>
#include "Math/Vector3D.h"
using V3 = ROOT::Math::XYZVector;

// Plane is defined by: (p - planePoint)·planeNormal = 0
// Ray is: p(t) = start + t * dir   (dir should be nonzero; unit is nice but not required)
// Returns intersection point. If no valid intersection, returns (NaN,NaN,NaN).
inline bool getIntersectionPoint(const V3 &planeNormal,
                                 const V3 &planePoint,
                                 const V3 &start_point,
                                 const V3 &direction,
                                 double &tOut,
                                 V3 &pOut,
                                 double eps = 1e-12)
{

    tOut = 0;

    double denom = planeNormal.Dot(direction);
    if (std::fabs(denom) < eps)
        return false;

    double t = planeNormal.Dot(planePoint - start_point) / denom;
    // printf ("t=%f\t",t);
    // If you only want forward intersections along the ray:
    if (t < 0)
        return false;
    tOut = t;
    pOut = start_point + t * direction;
    return true;
}

static inline V3 SampleIsotropic(TRandom3 &rng)
{
    double u = rng.Uniform(-1.0, 1.0);
    double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    double s = std::sqrt(std::max(0.0, 1.0 - u * u));
    return V3(s * std::cos(phi), s * std::sin(phi), u);
}

static inline V3 SampleInScint(TRandom3 &rng, double L, double W, double T, double wedgeLen, double wedgeTipW)
{
    double x = rng.Uniform(wedgeLen, L - wedgeLen);
    double y = rng.Uniform(-W / 2, +W / 2);
    double z = rng.Uniform(-T / 2, +T / 2);
    return V3(x, y, z);
}

void Transport_new()
{
    double L = 60.0;
    double W = 30.0;
    double T = 1.0;

    double wedgeLen = 10; //10.0;
    double wedgeTipW = 1;// 5.0;

    double zTop = 0.5 * T;
    double yTip = 0.5 * wedgeTipW;
    double yMid = 0.5 * W;
    // Use V3 everywhere
    std::vector<V3> pts;
    std::vector<V3> normals;
    std::vector<V3> pointPlane;

    V3 p_block_xNeg_yNeg_top(wedgeLen, -yMid, zTop);
    V3 p_block_xPos_yNeg_top(L - wedgeLen, -yMid, zTop);
    V3 p_block_xPos_yPos_top(L - wedgeLen, +yMid, zTop);
    V3 p_block_xNeg_yPos_top(wedgeLen, +yMid, zTop);
    V3 p_block_xNeg_yNeg_bottom(wedgeLen, -yMid, -zTop);
    V3 p_block_xPos_yNeg_bottom(L - wedgeLen, -yMid, -zTop);
    V3 p_block_xPos_yPos_bottom(L - wedgeLen, +yMid, -zTop);
    V3 p_block_xNeg_yPos_bottom(wedgeLen, +yMid, -zTop);
    V3 p_wedge_xPos_yNeg_bottom(L, -yTip, -zTop);
    V3 p_wedge_xPos_yPos_bottom(L, +yTip, -zTop);
    V3 p_wedge_xPos_yNeg_top(L, -yTip, zTop);
    V3 p_wedge_xPos_yPos_top(L, +yTip, zTop);
    V3 p_wedge_xNeg_yNeg_top(0, -yTip, zTop);
    V3 p_wedge_xNeg_yPos_top(0, +yTip, zTop);
    V3 p_wedge_xNeg_yNeg_bottom(0, -yTip, -zTop);
    V3 p_wedge_xNeg_yPos_bottom(0, +yTip, -zTop);


    V3 norm_top = -findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_block_xNeg_yNeg_top, p_wedge_xPos_yNeg_top);
    V3 norm_bottom = -findNormalFrom3Points(p_wedge_xPos_yNeg_bottom, p_block_xNeg_yNeg_bottom, p_wedge_xNeg_yNeg_bottom);
    V3 norm_yPos = -findNormalFrom3Points(p_block_xNeg_yPos_top, p_block_xPos_yPos_top, p_block_xNeg_yPos_bottom);
    V3 norm_yNeg = -findNormalFrom3Points(p_block_xNeg_yNeg_top, p_block_xNeg_yNeg_bottom, p_block_xPos_yNeg_bottom);
    V3 norm_xPos = -findNormalFrom3Points(p_wedge_xPos_yNeg_top, p_wedge_xPos_yNeg_bottom, p_wedge_xPos_yPos_bottom);
    V3 norm_xNeg = +findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_wedge_xNeg_yNeg_bottom, p_wedge_xNeg_yPos_bottom);
    V3 norm_wedge_ypos_xpos = +findNormalFrom3Points(p_wedge_xPos_yPos_top, p_block_xPos_yPos_top, p_block_xPos_yPos_bottom);
    V3 norm_wedge_yneg_xpos = -findNormalFrom3Points(p_wedge_xPos_yNeg_top, p_block_xPos_yNeg_top, p_block_xPos_yNeg_bottom);
    V3 norm_wedge_ypos_xneg = -findNormalFrom3Points(p_wedge_xNeg_yPos_top, p_block_xNeg_yPos_top, p_wedge_xNeg_yPos_bottom);
    V3 norm_wedge_yneg_xneg = +findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_block_xNeg_yNeg_top, p_wedge_xNeg_yNeg_bottom);
    // Plane 0: top
    normals.push_back(norm_top);
    pointPlane.push_back( p_block_xNeg_yNeg_top);
    // plane 1: bottom
    normals.push_back(norm_bottom);
    pointPlane.push_back( p_block_xNeg_yNeg_bottom);
    // plane 2: y positive
    normals.push_back(norm_yPos);
    pointPlane.push_back(p_block_xNeg_yPos_top);
    // plane 3: y negative
    normals.push_back(norm_yNeg);
    pointPlane.push_back(p_block_xPos_yNeg_top);
    // plane 4: x positive
    normals.push_back(norm_xPos);
    pointPlane.push_back(V3(L+2*wedgeLen, 0, 0));
    // plane 5 : x negative
    normals.push_back(norm_xNeg);
    pointPlane.push_back(V3(0,0,0));

    if (wedgeLen > 0 && wedgeTipW > 0)
    {
    //  plane 6: wedge y positive x positive
    normals.push_back(norm_wedge_ypos_xpos);
    pointPlane.push_back( p_wedge_xPos_yPos_top  );
    // plane 7: wedge y negative x positive
    normals.push_back(norm_wedge_yneg_xpos);
    pointPlane.push_back(  p_wedge_xPos_yNeg_top );
    // plane 8: wedge y positive x negative
    normals.push_back(norm_wedge_ypos_xneg);
    pointPlane.push_back(p_wedge_xNeg_yPos_top  );
    // plane 9: wedge y negative x negative
    normals.push_back(norm_wedge_yneg_xneg);
    pointPlane.push_back(p_wedge_xNeg_yNeg_top);
    }
    // Collect points for plotting
     for (int i = 0; i < (int)normals.size(); ++i)
     {
         const V3 &p = pointPlane[i];
         const V3 &n = normals[i];
         const V3 &offset = pointPlane[i];;// p + 0.1 * n;
         cout<<"pointPlane "<<i<<" at ("<<p.x()<<", "<<p.y()<<", "<<p.z()<<") with normal ("<<n.x()<<", "<<n.y()<<", "<<n.z()<<"):: ";
         cout << "i=" << i << "   is inside=" << isInside(offset, L, W, T, wedgeLen, wedgeTipW) << "  p=" << offset << "\n";
     }

    TRandom3 rng(0);
    int faceCounter[10] = {0};
     TGraph2D *ghits = new TGraph2D();
     ghits->SetMarkerColor(kRed);
     ghits->SetMarkerStyle(20);
    for (int j = 0; j < 1000 ; ++j)
    {
        V3 p0 = SampleInScint(rng, L, W, T, wedgeLen, wedgeTipW);
        V3 dir = SampleIsotropic(rng);
        printf("Starting point (%f, %f, %f), direction=(%f, %f, %f)\n", p0.x(), p0.y(), p0.z(), dir.x(), dir.y(), dir.z());

        int n = 0;
        int pre = 0;
        V3 prehitpoint;
        for (int i = 0; i < (int)normals.size(); ++i)
        { // cout<<"\n"<<i<<".\t point plabe="<<pointPlane[i]<<"   normal="<<normals[i]<<"\t";
            double t;
            V3 hitPoint;

            bool intersection = getIntersectionPoint(normals[i], pointPlane[i], p0, dir, t, hitPoint);

            if (intersection)
            {
               // if ((fabs(hitPoint.x()-L-2*wedgeLen)<1e-9 || hitPoint.x()>L+2*wedgeLen+1e-9 ) && (i==6 || i==7)) continue;
                cout<<hitPoint.x()<<"\t"<<hitPoint.y()<<"\t"<<hitPoint.z()<<"\n";
                int isin = isInside(hitPoint, L, W, T, wedgeLen, wedgeTipW);
                if (isin > 0)
                {
                    if (n == 0) {
                        pre = i;
                        prehitpoint = hitPoint;
                    }
                    if (n==1)  {
                        faceCounter[pre]++;
                         ghits->SetPoint(ghits->GetN(), prehitpoint.x(), prehitpoint.y(), prehitpoint.z());
                         ghits->SetPoint(ghits->GetN(), hitPoint.x(), hitPoint.y(), hitPoint.z());
                    }
                    if (n>0)
                    {
                            cout << "### Warning: double hit on same face " << i << " and " << pre << "\n";
                            faceCounter[i]++;
                    }
                    if (n > 0)
                        cout << "* ";
                    n++;
                    cout << n << ".   is inside=" << isin << " face =" << i << "\t t=" << t << "  hit=("
                         << hitPoint.x() << "," << hitPoint.y() << "," << hitPoint.z() << ")\n";
                }
            }

            // std::cout << "Intersection with plane " << i << " at (" << intersection.x() << ", " << intersection.y() << ", " << intersection.z() << ")";
        }
    }

    int tot = 0;
    for (int i = 0; i < 10; ++i)
    {
        cout << "Face type " << i << " was hit " << faceCounter[i] << " times.\n";
        tot += faceCounter[i];
    }
    cout << "Total hits=" << tot << "\n";

    pts.reserve(8);

    pts.push_back(V3(0.0, -yTip, zTop));
    pts.push_back(V3(wedgeLen, -yMid, zTop));
    pts.push_back(V3(L - wedgeLen, -yMid, zTop));
    pts.push_back(V3(L, -yTip, zTop));

    pts.push_back(V3(L, +yTip, zTop));
    pts.push_back(V3(L - wedgeLen, +yMid, zTop));
    pts.push_back(V3(wedgeLen, +yMid, zTop));
    pts.push_back(V3(0.0, +yTip, zTop));

    /* // normal sanity
     // V3 norm_top = findNormalFrom3Points(pts[0], pts[1], pts[3]);
     V3 testPoint1(p_wedge_xNeg_yNeg_top + 0.1 * norm_top);
     std::cout << "isInside(testPoint1)=" << isInside(testPoint1, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint2(p_wedge_xNeg_yNeg_bottom + 0.1 * norm_bottom);
     std::cout << "isInside(testPoint2)=" << isInside(testPoint2, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint3(p_wedge_xNeg_yPos_top + 0.1 * norm_yPos);
     std::cout << "isInside(testPoint3)=" << isInside(testPoint3, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint4(p_wedge_xNeg_yNeg_top + 0.1 * norm_yNeg);
     std::cout << "isInside(testPoint4)=" << isInside(testPoint4, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint5(p_wedge_xPos_yNeg_top + 0.1 * norm_xPos);
     std::cout << "isInside(testPoint5)=" << isInside(testPoint5, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint6(p_wedge_xNeg_yNeg_top + 0.1 * norm_xNeg);
     std::cout << "isInside(testPoint6)=" << isInside(testPoint6, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint7(p_wedge_xPos_yPos_top + 0.1 * norm_wedge_ypos_xpos);
     std::cout << "isInside(testPoint7)=" << isInside(testPoint7, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint8(p_wedge_xPos_yNeg_top + 0.1 * norm_wedge_yneg_xpos);
     std::cout << "isInside(testPoint8)=" << isInside(testPoint8, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint9(p_wedge_xNeg_yPos_top + 0.1 * norm_wedge_ypos_xneg);
     std::cout << "isInside(testPoint9)=" << isInside(testPoint9, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint10(p_wedge_xNeg_yNeg_top + 0.1 * norm_wedge_yneg_xneg);
     std::cout << "isInside(testPoint10)=" << isInside(testPoint10, L, W, T, wedgeLen, wedgeTipW) << "\n";
     V3 testPoint11((L + wedgeLen * 2) / 2, 0, 0);
     std::cout << "isInside(testPoint11)=" << isInside(testPoint11, L, W, T, wedgeLen, wedgeTipW) << "\n";*/
    // plot x-y
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TGraph2D *gpoints = new TGraph2D();

// ---------- TOP loop (8 points + close) ----------
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_top.x(),  p_wedge_xNeg_yNeg_top.y(),  p_wedge_xNeg_yNeg_top.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yNeg_top.x(),  p_block_xNeg_yNeg_top.y(),  p_block_xNeg_yNeg_top.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yNeg_top.x(),  p_block_xPos_yNeg_top.y(),  p_block_xPos_yNeg_top.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yNeg_top.x(),  p_wedge_xPos_yNeg_top.y(),  p_wedge_xPos_yNeg_top.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yPos_top.x(),  p_wedge_xPos_yPos_top.y(),  p_wedge_xPos_yPos_top.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yPos_top.x(),  p_block_xPos_yPos_top.y(),  p_block_xPos_yPos_top.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yPos_top.x(),  p_block_xNeg_yPos_top.y(),  p_block_xNeg_yPos_top.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yPos_top.x(),  p_wedge_xNeg_yPos_top.y(),  p_wedge_xNeg_yPos_top.z());
// close top loop
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_top.x(),  p_wedge_xNeg_yNeg_top.y(),  p_wedge_xNeg_yNeg_top.z());


// ---------- BOTTOM loop (8 points + close) ----------
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_bottom.x(),  p_wedge_xNeg_yNeg_bottom.y(),  p_wedge_xNeg_yNeg_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yNeg_bottom.x(),  p_block_xNeg_yNeg_bottom.y(),  p_block_xNeg_yNeg_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yNeg_bottom.x(),  p_block_xPos_yNeg_bottom.y(),  p_block_xPos_yNeg_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yNeg_bottom.x(),  p_wedge_xPos_yNeg_bottom.y(),  p_wedge_xPos_yNeg_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yPos_bottom.x(),  p_wedge_xPos_yPos_bottom.y(),  p_wedge_xPos_yPos_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yPos_bottom.x(),  p_block_xPos_yPos_bottom.y(),  p_block_xPos_yPos_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yPos_bottom.x(),  p_block_xNeg_yPos_bottom.y(),  p_block_xNeg_yPos_bottom.z());
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yPos_bottom.x(),  p_wedge_xNeg_yPos_bottom.y(),  p_wedge_xNeg_yPos_bottom.z());
// close bottom loop
gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_bottom.x(),  p_wedge_xNeg_yNeg_bottom.y(),  p_wedge_xNeg_yNeg_bottom.z());

/* for (int i = 0; i < (int)pts.size(); ++i)
    {
        const V3 &p = pts[i];
        // std::cout << "Corner " << i << " at (" << p.x() << ", " << p.y() << ", " << p.z() << ")\n";
        gpoints->SetPoint(i, p.x(), p.y(), p.z());
    }
        */

  //  gpoints->SetPoint((int)pts.size(), pts[0].x(), pts[0].y());

    gpoints->SetTitle("Top face outline; x [cm]; y [cm]");
    gpoints->SetMarkerStyle(20);
    gpoints->SetLineWidth(2);
    gpoints->Draw("AP LINE");
    ghits->Draw("P same");
    c->Update();
}
