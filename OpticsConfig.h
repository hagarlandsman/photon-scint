#ifndef OPTICSCONFIG_H
#define OPTICSCONFIG_H
#include <cmath>
#include "Math/Vector3D.h"
#include "TGraph2D.h"
using V3 = ROOT::Math::XYZVector;

#include <vector>


// -----------------------------
// Config and per-photon output
// -----------------------------
struct OpticsConfig
{
    double nScint = 1.58;
    double nOut = 1.0;

    double absLen = 300.0;
    double Rwrap = 0.95;
    int maxSteps = 2000;

    double rPMT = 2.5; // cm
    double epsCouple = 0.90;
    double pde = 0.20;

    double L = 90 ;// cm
    double W = 30 ; // cm
    double T = 1; //cm
    int wrap = 1; // 1=PTFE, 2=Mylar


    // coupling position dependence
    double eps0 = 0.00;     // Kills all indirect hits when =0
    double lambdaC = 120.0; // irrelvant when eps=0

    // toggles
    bool savePath = false;
    bool mirrorx = true;
    bool mirrory = true;
    double dScan = 2.0; // cm
    // wedge parameters
    bool useWedge = false;
    double wedgeLen = 20.0; // cm
    double wedgeTipW = 5.0; // cm (width at end plane)
};

struct PhotonResult
{
    // stamped inputs
    double L = 0, W = 0, T = 0;
    int wrap = 1; // 1=PTFE, 2=Mylar
    int site_number = 0;

    // site + end point
    double x0 = 0, y0 = 0, z0 = 0;
    double xf = 0, yf = 0, zf = 0;

    // outcome
    int endPlane = -1; // 0..5 box, 6/7 wedge walls
    int pmt_side = 0;  // 1=x=0, 2=x=L
    int inPMT = 0;
    int detected = 0;
    int absorbed = 0;
    int escaped = 0;
    int reachedEnd = 0;

    int nBounces = 0;
    double path = 0;
    double epsRaysPMT = 0;

    // optional path
    std::vector<float> xPath, yPath, zPath;
};


using V3 = ROOT::Math::XYZVector;
inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double eps = 1e-6);

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double wedgeLen, double wedgeTipW,
                    double eps = 1e-6);

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double eps)
{
    return isInside(p, L, W, T, 0.0, 0.0, eps);
}

inline int isInside(const V3 &p, const OpticsConfig &cfg) {
                    double L=cfg.L;
                    double W=cfg.W;
                    double T=cfg.T;
                    double wedgeLen=cfg.useWedge ? cfg.wedgeLen : 0.0;
                    double wedgeTipW=cfg.useWedge ? cfg.wedgeTipW : 0.0;
                    double eps=1e-6;

    return isInside(p, L, W, T, wedgeLen, wedgeTipW, eps);
}

inline int isInside(const V3 &p,
                    double L, double W, double T,
                    double wedgeLen, double wedgeTipW,
                    double eps)
{

    // x,y,z check
    if (fabs(p.x()) > L * 0.5 + wedgeLen + eps) {
       // cout<<" isInside: x="<<p.x()<<" out of bounds ("<<-L*0.5 - wedgeLen - eps<<","<<L*0.5 + wedgeLen + eps<<") \n";
        return 0;
    }
    if (fabs(p.y()) > W * 0.5 + eps) {
     //cout<<" isInside: y="<<p.y()<<" out of bounds ("<<-W*0.5 - eps<<","<<W*0.5 + eps<<") \n";
        return 0;
    }
    if (fabs(p.z()) > T * 0.5 + eps)
    {
   //  cout<<" isInside: z="<<p.z()<<" out of bounds ("<<-T*0.5 - eps<<","<<T*0.5 + eps<<") \n";
        return 0;
    }
    // Main block - if in central part return 1
    if (fabs(p.x()) <= L * 0.5 + eps && fabs(p.y()) <= 0.5 * W + eps && fabs(p.z()) <= 0.5 * T + eps)
        return 1;

    if (fabs(p.x()) <= L * 0.5 + wedgeLen + eps)
    {
        double y2 = 0.5 * W;
        double y1 = 0.5 * wedgeTipW;
        double x2 = L * 0.5;
        double x1 = L * 0.5 + wedgeLen;
        double slope = (y2 - y1) / (x2 - x1); // negative
        double yMax = y1 + slope * (fabs(p.x()) - x1);
        if (fabs(p.y()) <= yMax + eps)
            return 2;
        else
        {
          //  cout<<"eps="<<eps<<"\n";
           // cout<<" isInside: x="<<p.x()<<" within wedge x-range but y=fabs("<<p.y()<<") exceeds wedge yMax="<<yMax<<" Delta="<<fabs(p.y()) - yMax<<"\n";
            return 0;
        }
    }
    return 0; // not in main block, but x is within main block x-range, so must be outside
}
// Normal from 3 points: (p1-p0) x (p2-p0), normalized
static V3 findNormalFrom3Points(const V3 &p0, const V3 &p1, const V3 &p2)
{
    V3 n = (p1 - p0).Cross(p2 - p0);
    double mag = std::sqrt(n.Dot(n));
    if (mag > 0)
        n = (1.0 / mag) * n;
    double epsilon = 1e-6;
    if (std::fabs(n.x()) < epsilon)
        n.SetX(0.0);
    if (std::fabs(n.y()) < epsilon)
        n.SetY(0.0);
    if (std::fabs(n.z()) < epsilon)
        n.SetZ(0.0);
    // std::cout << "n = (" << n.x() << ", " << n.y() << ", " << n.z() << ")\n";
    return n;
}
void makeFaces(const OpticsConfig& cfg,
               TGraph2D*& gpoints,
               std::vector<V3>& normals,
               std::vector<V3>& pointPlane)
{
    // reset outputs
    if (gpoints) { delete gpoints; gpoints = nullptr; }
    normals.clear();
    pointPlane.clear();

    double L = cfg.L;
    double W = cfg.W;
    double T = cfg.T;
    double wedgeLen  = cfg.useWedge ? cfg.wedgeLen  : 0.0;
    double wedgeTipW = cfg.useWedge ? cfg.wedgeTipW : 0.0;

    double yMid = W * 0.5;
    double yTip = wedgeTipW * 0.5;
    double zTop = T * 0.5;

    // (left-right centered in x): your exact point definitions
    V3 p_wedge_xNeg_yNeg_top(-L * 0.5 - wedgeLen, -yTip, zTop);
    V3 p_block_xNeg_yNeg_top(-L * 0.5, -yMid, zTop);
    V3 p_block_xPos_yNeg_top(+L * 0.5, -yMid, zTop);
    V3 p_wedge_xPos_yNeg_top(+L * 0.5 + wedgeLen, -yTip, zTop);
    V3 p_wedge_xPos_yPos_top(+L * 0.5 + wedgeLen, +yTip, zTop);
    V3 p_block_xPos_yPos_top(+L * 0.5, yMid, zTop);
    V3 p_block_xNeg_yPos_top(-L * 0.5, +yMid, zTop);
    V3 p_wedge_xNeg_yPos_top(-L * 0.5 - wedgeLen, +yTip, zTop);

    V3 p_wedge_xNeg_yNeg_bottom(-L * 0.5 - wedgeLen, -yTip, -zTop);
    V3 p_block_xNeg_yNeg_bottom(-L * 0.5, -yMid, -zTop);
    V3 p_block_xPos_yNeg_bottom(+L * 0.5, -yMid, -zTop);
    V3 p_wedge_xPos_yNeg_bottom(+L * 0.5 + wedgeLen, -yTip, -zTop);
    V3 p_wedge_xPos_yPos_bottom(+L * 0.5 + wedgeLen, +yTip, -zTop);
    V3 p_block_xPos_yPos_bottom(+L * 0.5, yMid, -zTop);
    V3 p_block_xNeg_yPos_bottom(-L * 0.5, +yMid, -zTop);
    V3 p_wedge_xNeg_yPos_bottom(-L * 0.5 - wedgeLen, +yTip, -zTop);

    // graph output
    gpoints = new TGraph2D();

    // ---------- TOP loop (8 points + close) ----------
    gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_top.x(),  p_wedge_xNeg_yNeg_top.y(),  p_wedge_xNeg_yNeg_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yNeg_top.x(),  p_block_xNeg_yNeg_top.y(),  p_block_xNeg_yNeg_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yNeg_top.x(),  p_block_xPos_yNeg_top.y(),  p_block_xPos_yNeg_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yNeg_top.x(),  p_wedge_xPos_yNeg_top.y(),  p_wedge_xPos_yNeg_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_wedge_xPos_yPos_top.x(),  p_wedge_xPos_yPos_top.y(),  p_wedge_xPos_yPos_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_block_xPos_yPos_top.x(),  p_block_xPos_yPos_top.y(),  p_block_xPos_yPos_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_block_xNeg_yPos_top.x(),  p_block_xNeg_yPos_top.y(),  p_block_xNeg_yPos_top.z());
    gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yPos_top.x(),  p_wedge_xNeg_yPos_top.y(),  p_wedge_xNeg_yPos_top.z());
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
    gpoints->SetPoint(gpoints->GetN(), p_wedge_xNeg_yNeg_bottom.x(),  p_wedge_xNeg_yNeg_bottom.y(),  p_wedge_xNeg_yNeg_bottom.z());

    // normals (unchanged)
    V3 norm_top              = +findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_block_xNeg_yNeg_top, p_wedge_xPos_yNeg_top);
    V3 norm_bottom           = +findNormalFrom3Points(p_wedge_xPos_yNeg_bottom, p_block_xNeg_yNeg_bottom, p_wedge_xNeg_yNeg_bottom);
    V3 norm_yPos             = +findNormalFrom3Points(p_block_xNeg_yPos_top, p_block_xPos_yPos_top, p_block_xNeg_yPos_bottom);
    V3 norm_yNeg             = +findNormalFrom3Points(p_block_xNeg_yNeg_top, p_block_xNeg_yNeg_bottom, p_block_xPos_yNeg_bottom);
    V3 norm_xPos             = -findNormalFrom3Points(p_wedge_xPos_yNeg_top, p_wedge_xPos_yNeg_bottom, p_wedge_xPos_yPos_bottom);
    V3 norm_xNeg             = +findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_wedge_xNeg_yNeg_bottom, p_wedge_xNeg_yPos_bottom);
    V3 norm_wedge_ypos_xpos  = +findNormalFrom3Points(p_wedge_xPos_yPos_top, p_block_xPos_yPos_top, p_block_xPos_yPos_bottom);
    V3 norm_wedge_yneg_xpos  = -findNormalFrom3Points(p_wedge_xPos_yNeg_top, p_block_xPos_yNeg_top, p_block_xPos_yNeg_bottom);
    V3 norm_wedge_ypos_xneg  = -findNormalFrom3Points(p_wedge_xNeg_yPos_top, p_block_xNeg_yPos_top, p_wedge_xNeg_yPos_bottom);
    V3 norm_wedge_yneg_xneg  = +findNormalFrom3Points(p_wedge_xNeg_yNeg_top, p_block_xNeg_yNeg_top, p_wedge_xNeg_yNeg_bottom);

    // planes 0..5 (unchanged)
//    normals.push_back(norm_xPos);
    normals.push_back(V3(-1,0,0));
    pointPlane.push_back(V3(L * 0.5 + wedgeLen, 0, 0));

//    normals.push_back(norm_xNeg);
    normals.push_back(V3(1,0,0));
    pointPlane.push_back(V3(-L * 0.5 - wedgeLen, 0, 0));


   // normals.push_back(norm_top);
    normals.push_back(V3(0,0,-1));
   pointPlane.push_back(p_block_xNeg_yNeg_top);

//    normals.push_back(norm_bottom);
    normals.push_back(V3(0,0,1));
    pointPlane.push_back(p_block_xNeg_yNeg_bottom);

  //  normals.push_back(norm_yPos);
    normals.push_back(V3(0,-1,0));
    pointPlane.push_back(p_block_xNeg_yPos_top);

  //  normals.push_back(norm_yNeg);
    normals.push_back(V3(0,1,0));
    pointPlane.push_back(p_block_xPos_yNeg_top);

    // wedges 6..9 (unchanged gating)
    if (wedgeLen > 0 && wedgeTipW > 0)
    {
        normals.push_back(norm_wedge_ypos_xpos);
        pointPlane.push_back(p_wedge_xPos_yPos_top);

        normals.push_back(norm_wedge_yneg_xpos);
        pointPlane.push_back(p_wedge_xPos_yNeg_top);

        normals.push_back(norm_wedge_ypos_xneg);
        pointPlane.push_back(p_wedge_xNeg_yPos_top);

        normals.push_back(norm_wedge_yneg_xneg);
        pointPlane.push_back(p_wedge_xNeg_yNeg_top);
    }
  /*   for (size_t i = 0; i < normals.size(); i++)
    {
        const V3 &p = pointPlane[i];
        const V3 &n = normals[i];
        cout << "pointPlane " << i << " at (" << p.x() << ", " << p.y() << ", " << p.z() << ") with normal (" << n.x() << ", " << n.y() << ", " << n.z() << ")\n";
    }
*/
}



#endif
