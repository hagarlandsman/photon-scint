/*
.L TreeWriter.cxx+
.L Geometry.cxx+
.L Transport.cxx+
.L RunToy.cxx+

.L Geometry.cxx+
.L Transport.cxx+
.L TreeWriter.cxx+
.L ToyOpticsWedge.C+
    OpticsConfig cfg;
    cfg.savePath = true;

    // geometry (also passed explicitly below, matching your current PropagateOnePhoton signature)
    cfg.L = 90.0;
    cfg.W = 30.0;
    cfg.T = 1.0;
    cfg.wrap = 1; // 1=PTFE, 2=Mylar

    // optics
    cfg.nScint = 1.58;
    cfg.nOut = 1.0;
    cfg.absLen = 300.0;
    cfg.Rwrap = 0.95;
    cfg.maxSteps = 2000;

    // PMT / coupling
    cfg.rPMT = 2.5;
    cfg.epsCouple = 0.90;
    cfg.pde = 0.20;
    cfg.eps0 = 0.0;
    cfg.lambdaC = 120.0;

    // wedges
    cfg.useWedge = true;
    cfg.wedgeLen = 20.0;
    cfg.wedgeTipW = 5.0;

    const char *outFile = "onePoint_100.root";
   RunManySites(10,10, cfg, outFile)
    DrawEvent4ViewFromTree("onePoint_100.root",1,true)
    DrawEventFromTree("onePoint_100.root",1,true,false)
*/
#include "TRandom3.h"
#include <iostream>
#include "Geometry.h"
#include "Vec3.h"
#include "OpticsConfig.h"
#include "TreeWriter.h"
#include "Transport.h"

void RunOnePoint100AndDraw()

{
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
    gROOT->ProcessLine(".L RunToy.cxx+");
    gROOT->ProcessLine(".L DrawHelper.cxx+");





    OpticsConfig cfg;
    cfg.savePath = true;

    // geometry (also passed explicitly below, matching your current PropagateOnePhoton signature)
    cfg.L = 90.0;
    cfg.W = 30.0;
    cfg.T = 1.0;
    cfg.wrap = 1; // 1=PTFE, 2=Mylar

    // optics
    cfg.nScint = 1.58;
    cfg.nOut = 1.0;
    cfg.absLen = 300.0;
    cfg.Rwrap = 0.95;
    cfg.maxSteps = 2000;

    // PMT / coupling
    cfg.rPMT = 2.5;
    cfg.epsCouple = 0.90;
    cfg.pde = 0.20;
    cfg.eps0 = 0.0;
    cfg.lambdaC = 120.0;

    // wedges
    cfg.useWedge = true;
    cfg.wedgeLen = 20.0;
    cfg.wedgeTipW = 5.0;

    const char *outFile = "onePoint_100.root";

    // fixed emission point (pick anything inside the active volume)
    Vec3 site(cfg.L * 0.5, 0.0, 0.0);

    TRandom3 rng(0);
    TreeWriter wr(outFile, cfg);

    const int N = 100;
    for (int i = 0; i < N; i++)
    {
        PhotonResult res = PropagateOnePhoton(
            rng,
            site,
            0,               // site_number
            cfg);

        wr.Fill(res);
    }

    wr.Close();

    // Draw the 1st event
  /*  DrawEventFromTree("wedge.root", 0, true, false);
DrawEventSplitViewFromTree("wedge.root", 1, true, 20.0, 5.0, 0.10);
DrawEvent4ViewFromTree("wedge.root",   1, true,  20.0, 5.0, 0.12);
// Zoom in tighter:
DrawEventFromTree("wedge.root", 0, true, 20.0, 5.0, true, 0.05);
*/

    //DrawEventFromTree(outFile, 0, true,  true);
}







