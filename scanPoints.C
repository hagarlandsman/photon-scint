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
#include "DrawHelper.h"
void scanPoints()

{
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
    //  gROOT->ProcessLine(".L RunToy.cxx+");
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

    double x0, x1, y0, y1, z0, z1;
    if (cfg.useWedge)
    {
        x0 = cfg.wedgeLen;
        x1 = cfg.L - cfg.wedgeLen;
    }
    else if (!cfg.useWedge)
    {
        x0 = 0;
        x1 = cfg.L;
    }
    y0 = -cfg.W * 0.5;
    y1 = cfg.W * 0.5;
    z0 = -cfg.T * 0.5;
    ;
    z1 = cfg.T * 0.5;
    ;
    double epsilon = 1e-5; // 1e-6;
    // Vec3 site((cfg.L) * 0.5, 0, 0); // center
    //  Vec3 site(x0 + epsilon ,y0 + epsilon,z0 + epsilon); // corner 1
   // Vec3 site(x1 - epsilon, y1 - epsilon, z1 - epsilon); // corner 2
    // printf ("%f %f %f\n",site.x,site.y,site.z);

    int NstepsX = 3;
    int NstepsY = 3;

    for (int ix = 0; ix < NstepsX; ix++)
    {
        double x = x0 + (x1 - x0) * ix / (NstepsX - 1);
        for (int iy = 0; iy < NstepsY; iy++)
        {
            double y = y0 + (y1 - y0) * iy / (NstepsY - 1);
            TString outFile = Form("out_%03d_%03d.root", ix, iy); //    const char *outFile = "onePoint_100.root";
            TRandom3 rng(0);
            TreeWriter wr(outFile, cfg);
            Vec3 site(x, y, 0); // corner 2
            printf ("%d %d %f %f \n",ix,iy,x,y);
            const int N = 100;
            for (int i = 0; i < N; i++)
            {
                PhotonResult res = PropagateOnePhoton(
                    rng,
                    site,
                    0, // site_number
                    cfg);

                wr.Fill(res);
            }

            wr.Close();
            DrawEventSplitViewFromTree(outFile, 12, true, 0.12);
         //   DrawEvent4ViewFromTree("onePoint_100.root", 12, true, 0.12);
        }
    }
            // Draw the 1st event
            /*  DrawEventFromTree("wedge.root", 0, true, false);
          DrawEventSplitViewFromTree("wedge.root", 1, true, 20.0, 5.0, 0.10);
          DrawEvent4ViewFromTree("wedge.root",   1, true,  20.0, 5.0, 0.12);
          // Zoom in tighter:
          DrawEventFromTree("wedge.root", 0, true, 20.0, 5.0, true, 0.05);
          */

            // DrawEventFromTree(outFile, 0, true,  true);
        }
