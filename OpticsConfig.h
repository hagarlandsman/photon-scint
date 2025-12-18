#ifndef OPTICSCONFIG_H
#define OPTICSCONFIG_H

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



#endif
