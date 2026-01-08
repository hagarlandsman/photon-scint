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
//
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPolyLine3D.h"
#include "TGraph2D.h"
#include "TPolyLine.h"
#include "TStyle.h"


//
#include "Math/Vector3D.h"


using V3 = ROOT::Math::XYZVector;

static inline V3 SampleInScint(TRandom3 &rng, double L, double W, double T, double wedgeLen, double wedgeTipW)
{

    double x = rng.Uniform(-L * 0.5, +L * 0.5);
    double y = rng.Uniform(-W * 0.5, +W * 0.5);
    double z = rng.Uniform(-T * 0.5, +T * 0.5);
    return V3(x, y, z);
}
void RunOnePoint100AndDraw(double WW_=30.0, double LL_=90, double TT_=1, int withwedge_=0, int wrap=1)

{
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
  //  gROOT->ProcessLine(".L RunToy.cxx+");
    gROOT->ProcessLine(".L DrawHelper.cxx+");
//   DrawEvent4ViewFromTree("onePoint_100.root",   12, true,     0.12);

    OpticsConfig cfg;
    cfg.savePath = false;

    // geometry (also passed explicitly below, matching your current PropagateOnePhoton signature)
    cfg.L = LL_; //90.0;
    cfg.W = WW_ ; //30.0;
    cfg.T = TT_; //1.0;
    //cfg.T = 1.0;
    cfg.wrap = wrap; // 1=PTFE, 2=Mylar

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
    if (withwedge_==0) {
        cout<<"No wedge\n";
        cfg.useWedge = false;
        cfg.wedgeLen = 0;
        cfg.wedgeTipW = 0;
    } else
    if (withwedge_==1) {
        cout<<"With wedge\n";
        cfg.useWedge = true;
        cfg.wedgeLen = 20; // 20.0;
        cfg.wedgeTipW = 5; // 5.0;
    }
     if (withwedge_==2) {
        cout<<"With wedge\n";
        cfg.useWedge = true;
        cfg.wedgeLen = WW_ * 0.5; // 20.0;
        cfg.wedgeTipW = 5; // 5.0;
    }
    //cfg.useWedge = true;
    //cfg.wedgeLen = 25; // 20.0;
    //cfg.wedgeTipW = 5; // 5.0;

    const char *outFile = Form("randPoint_W%d_L%d_T%d_Wedge%d_Wrap%d.root", (int)WW_, (int)LL_, (int)TT_, withwedge_, wrap);

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

    TPolyLine3D* pl = BuildFrameGeometry3D(g2);
    pl->SetLineColor(kBlue);
    pl->SetLineWidth(2);
    TCanvas* cgeo = (TCanvas*)gROOT->FindObject("cgeo");
    if (!cgeo) cgeo = new TCanvas("cgeo","Scintillator Geometry",800,600);
    cgeo->cd();
    cgeo->Clear();
    g2->Draw("P0");       // or your 3D frame that creates axes

    pl->Draw("same");
    double xrange = -1.2*(cfg.L * 0.5 + cfg.wedgeLen);
    double yrange = 1.2*(cfg.W * 0.5);
    double zrange = 1.2*(cfg.T * 0.5);
    SetViewEqualXYZ(-xrange, xrange, -yrange, yrange, -zrange, zrange, true);

   // g2->Write("scintGeometry");
   cgeo->Modified();
    cgeo->Update();
    gSystem->ProcessEvents();
    const int N = 100 * cfg.L * cfg.W ; // 200000;
    cout<<"W="<<cfg.W<<" L="<<cfg.L<<" Running "<<N<<" photons...\n";
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


    wr.Close();
  // DrawEventSplitViewFromTree("onePoint_100.root", 12, true, 0.12);
  if (cfg.savePath) {
  DrawEvent4ViewFromTree("onePoint_100.root",   0,      0.12);
   DrawEventFromTree3D("onePoint_100.root",0);
  }
    const double cost_cm2 = 0.4; // Posswibly extended to cm^3 with cfg.T=1.0 cm
    const double cost_pmt_1inch = 1269.0;
    const double cost_pmt_2inch = 1890.0;

    double costPMT;
    if (cfg.rPMT <= 1.27)
        costPMT = cost_pmt_1inch;
    else if (cfg.rPMT <= 2.54)
        costPMT = cost_pmt_2inch;
    double cost = cfg.T*cfg.W*cfg.L*cost_cm2 + 2 * costPMT;

    double fDet = 1.*Ndet / Ntot;
    // in this simulation each scintillator get the same number of photons
    // A larger scintillator with same efficiency will detect more photons in total

    double fDet_total = fDet * (cfg.W * cfg.L);
    double cost_per_eff_area = cost / (fDet * cfg.W * cfg.L);

    cout<<"Done "<<N<<endl;
    cout<< "Nabs: "<<1.*Nabs / Ntot *100. <<"% \n";
    cout << "Nesc:" << 1.* Nesc / Ntot * 100.<<"% \n";
    cout << "NhPMT:" << 1.* NhPMT / Ntot * 100.<<"% \n";
    cout << "Ndet:" << 1.*Ndet / Ntot * 100.<<"% /pm "<< sqrt( 1.*Ndet*(Ntot-Ndet)/Ntot/Ntot/Ntot)*100.<<"%\n";

    cout << "Detected fraction per scintillator: " << fDet * 100.0 << "%\n";
    cout << "Detected fraction total: " << fDet_total * 100.0 << "%\n";
    cout<< "Cost estimate: $" << cost <<endl
        << "\t Area=" << (cfg.W * cfg.L) <<"at "<<cost_cm2<<"$/cm^2"<<endl
        << "\t 2 PMTs r="<<cfg.rPMT << "cm, at " << costPMT << "$ each"<<endl;
    cout<< "real area = "<< (cfg.W*cfg.L) << " cm^2\n";
    cout<< "effective area = "<< fDet * (cfg.W*cfg.L) << " cm^2\n";
    cout<< "cost per effective="<< cost_per_eff_area << "\n";
    cout<<" Geometry: L="<<cfg.L<<", W="<<cfg.W<<", T="<<cfg.T<<"\n";
    cout<<" Wedge: useWedge="<<cfg.useWedge<<", wedgeLen="<<cfg.wedgeLen<<", wedgeTipW="<<cfg.wedgeTipW<<"\n";
    cout<<" Optics: nScint="<<cfg.nScint<<", nOut="<<cfg.nOut<<", absLen="<<cfg.absLen<<", Rwrap="<<cfg.Rwrap << "\n";
    cout<<" PMT/coupling: rPMT="<<cfg.rPMT<<", epsCouple="<<cfg.epsCouple<<", pde="<<cfg.pde<<", eps0="<<cfg.eps0<<", lambdaC="<<cfg.lambdaC<<"\n";
    cout<<" Wrap: "<<((cfg.wrap==1)?"PTFE":"Mylar")<<"\n";
    printf ("Ntot,Ndet/Ntot,err,Nabs/Ntot,Nesc/Ntot,NhPMT/Ntot,cost,fDet_total,cost_per_eff_area,L,W,T,useWedge,wedgeLen,wedgeTipW,nScint,absLen,Rwrap\n");
    printf ("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        Ntot*1.,
        1.*Ndet / Ntot,
        sqrt( 1.*Ndet*(Ntot-Ndet)/Ntot/Ntot/Ntot),
        1.*Nabs / Ntot,
        1.*Nesc / Ntot,
        1.*NhPMT / Ntot,
        cost,
        fDet_total,
        cost_per_eff_area,
        cfg.L,
        cfg.W,
        cfg.T,
        cfg.useWedge?1.0:0.0,
        cfg.wedgeLen,
        cfg.wedgeTipW,
        cfg.nScint,
        cfg.absLen,
        cfg.Rwrap
    );
    TString fout_name="random_points.csv";
    FILE *f = fopen (fout_name.Data(),"a");
    fprintf(f,"Ntot,Ndet/Ntot,err,Nabs/Ntot,Nesc/Ntot,NhPMT/Ntot,cost,fDet_total,cost_per_eff_area,L,W,T,useWedge,wedgeLen,wedgeTipW,nScint,absLen,Rwrap,wrap,rPMT\n");
    fprintf (f,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f\n",
        Ntot*1.,
        1.*Ndet / Ntot,
        sqrt( 1.*Ndet*(Ntot-Ndet)/Ntot/Ntot/Ntot),
        1.*Nabs / Ntot,
        1.*Nesc / Ntot,
        1.*NhPMT / Ntot,
        cost,
        fDet_total,
        cost_per_eff_area,
        cfg.L,
        cfg.W,
        cfg.T,
        cfg.useWedge?1.0:0.0,
        cfg.wedgeLen,
        cfg.wedgeTipW,
        cfg.nScint,
        cfg.absLen,
        cfg.Rwrap,
        cfg.wrap,
        cfg.rPMT);
    fclose (f);
   // Draw the 1st event
  /*  DrawEventFromTree("wedge.root", 0, true, false);
DrawEventSplitViewFromTree("wedge.root", 1, true, 20.0, 5.0, 0.10);
DrawEvent4ViewFromTree("wedge.root",   1, true,  20.0, 5.0, 0.12);
// Zoom in tighter:
DrawEventFromTree("wedge.root", 0, true, 20.0, 5.0, true, 0.05);
*/

    //DrawEventFromTree(outFile, 0, true,  true);
}







