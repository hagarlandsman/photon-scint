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
InsideActiveVolume
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
//#include "Vec3.h"
#include "OpticsConfig.h"
#include "TreeWriter.h"
#include "Transport.h"
#include "DrawHelper.h"

#include "Math/Vector3D.h"

using V3 = ROOT::Math::XYZVector;

static inline V3 SampleInScint(TRandom3 &rng, double L, double W, double T, double wedgeLen, double wedgeTipW)
{

    double x = rng.Uniform(-L * 0.5, +L * 0.5);
    double y = rng.Uniform(-W * 0.5, +W * 0.5);
    double z = rng.Uniform(-T * 0.5, +T * 0.5);
    return V3(x, y, z);
}
void RunOnePoint100AndDraw()

{
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
  //  gROOT->ProcessLine(".L RunToy.cxx+");
    gROOT->ProcessLine(".L DrawHelper.cxx+");
//   DrawEvent4ViewFromTree("onePoint_100.root",   12, true,     0.12);

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
    cfg.maxSteps = 4000;

    // PMT / coupling
    cfg.rPMT = 1.27;
    cfg.epsCouple = 0.90;
    cfg.pde = 0.20;
    cfg.eps0 = 0.0;
    cfg.lambdaC = 120.0;

    // wedges
    //cfg.useWedge = true;
    //cfg.wedgeLen = 25; // 20.0;
    //cfg.wedgeTipW = 5; // 5.0;
    cfg.useWedge = false;
    cfg.wedgeLen = 0; // 20.0;
    cfg.wedgeTipW = 0; // 5.0;

    const char *outFile = "onePoint_100.root";

    // fixed emission point (pick anything inside the active volume)

    double x0,x1,yy0,yy1,z0,z1;
    x0 = 0;
    x1 = cfg.L;
    if (cfg.useWedge && cfg.wedgeLen>0 ) {
        x0 = cfg.wedgeLen;
        x1 = cfg.L - cfg.wedgeLen;
    }


     // center
    TRandom3 rng(0);

    V3 site(0,0,0);

    TreeWriter wr(outFile, cfg);

    TTree *tCfg = new TTree("tCfg","Run configuration");
    tCfg->Branch("cfg", &cfg);
    tCfg->Fill();
    tCfg->Write();
    TGraph2D *g2 = nullptr;
    std::vector<V3> normals, pointPlane;
    cout << "Making faces...\n";
    makeFaces(cfg, g2, normals, pointPlane);
    cout << "Making faces doe...\n";
    g2->Draw("AP LINE");
    g2->Write("scintGeometry");
    const int N = 1000000 ;
    int Nabs = 0;
    int Ntot = 0;
    int Nesc = 0;
    int NhPMT = 0;
    int Ndet = 0;
    for (int i = 0; i < N; i++)
    {
            site = SampleInScint(rng, cfg.L, cfg.W, cfg.T, cfg.wedgeLen, cfg.wedgeTipW);
      //  cout<<"site="<<site<<endl;
       // cout<<"----------------------\ni="<<i<<"\n";
        PhotonResult res = PropagateOnePhoton(
            rng,
            site,
            0,               // site_number
            cfg,normals, pointPlane);

             Ntot++;
            Nabs += res.absorbed ;
            Nesc += res.escaped ;
            NhPMT += res.inPMT ;
            Ndet += res.detected ;
        wr.Fill(res);
    }
    cout<<"Done "<<N<<endl;
    cout<< "Nabs: "<<1.*Nabs / Ntot *100. <<"% \n";
    cout << "Nesc:" << 1.* Nesc / Ntot * 100.<<"% \n";
    cout << "NhPMT:" << 1.* NhPMT / Ntot * 100.<<"% \n";
    cout << "Ndet:" << 1.*Ndet / Ntot * 100.<<"% /pm "<< sqrt( 1.*Ndet*(Ntot-Ndet)/Ntot/Ntot/Ntot)*100.<<"%\n";

    wr.Close();
  // DrawEventSplitViewFromTree("onePoint_100.root", 12, true, 0.12);
   DrawEvent4ViewFromTree("onePoint_100.root",   0,      0.12);
   DrawEventFromTree3D("onePoint_100.root",0);
   // Draw the 1st event
  /*  DrawEventFromTree("wedge.root", 0, true, false);
DrawEventSplitViewFromTree("wedge.root", 1, true, 20.0, 5.0, 0.10);
DrawEvent4ViewFromTree("wedge.root",   1, true,  20.0, 5.0, 0.12);
// Zoom in tighter:
DrawEventFromTree("wedge.root", 0, true, 20.0, 5.0, true, 0.05);
*/

    //DrawEventFromTree(outFile, 0, true,  true);
}







