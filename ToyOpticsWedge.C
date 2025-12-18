// ToyOpticsWedge.C
// ROOT macro for toy optical photon transport in a rectangular scintillator bar
// with optional end light-guide wedges (taper in y only, no taper in z).
//
// Compile in ROOT:
//   .L ToyOpticsWedge.C+
//
// Example run (writes a tree):
//   OpticsConfig cfg;
//   cfg.savePath = true;      // store xPath/yPath/zPath
//   cfg.useWedge = true;      // enable wedges at both ends
//   cfg.wedgeLen = 20.0;      // cm
//   cfg.wedgeTipW = 5.0;      // cm (width at x=0 and x=L)
//   RunManySites(50, 200, 120, 20, 1, 1, cfg, "toyOptics.root");
//
// Draw one event from the file (if path was saved):
//   DrawEventFromTree("toyOptics.root", 0, true);

#include "TRandom3.h"
#include "Geometry.h"
#include "Vec3.h"
#include "OpticsConfig.h"
#include <iostream>
#include "TreeWriter.h"
#include "Transport.h"




// -----------------------------
// Optics + geometry helpers
// -----------------------------
// Linear wedge half-width yMax(x)
// -----------------------------
// Propagate one photon from a given site position
// -----------------------------

// -----------------------------
// Tree writer
// -----------------------------

// -----------------------------
// Wrapper: many sites and many photons per site
// -----------------------------
void RunManySites(
    long long Nsites,
    long long NphotPerSite,
    const OpticsConfig &cfg,
    const char *outFile = "toyOptics.root")
{
    TRandom3 rng(0);
    TreeWriter wr(outFile, cfg);
    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;
    long long nDet = 0;
    for (long long j = 0; j < Nsites; j++)
    {

        Vec3 site;

        do
        {
            site = Vec3(
                rng.Uniform(0.0, L),
                rng.Uniform(-W / 2.0, W / 2.0),
                rng.Uniform(-T / 2.0, T / 2.0));
        } while (!InsideActiveVolume(site.x, site.y, site.z,  cfg));

        for (long long i = 0; i < NphotPerSite; i++)
        {
            PhotonResult res = PropagateOnePhoton(rng, site, (int)j,  cfg);
            if (res.detected)
                nDet++;
            wr.Fill(res);
        }
    }

    wr.Close();

    double frac = (Nsites > 0 && NphotPerSite > 0) ? double(nDet) / double(Nsites * NphotPerSite) : 0.0;
    std::cout << "RunManySites: detected fraction = " << frac
              << " (Nsites=" << Nsites << ", NphotPerSite=" << NphotPerSite << ")\n";
}

void ScanBoard(   // HYL fix this - need to write with two loops
    long long Nsites,
    long long NphotPerSite,
    const OpticsConfig &cfg,
    const char *outFile = "toyOptics.root")
{
    TRandom3 rng(0);
    TreeWriter wr(outFile, cfg);
    double L = cfg.L;
    double T = cfg.T;
    double W = cfg.W;
    int wrap = cfg.wrap;
    long long nDet = 0;
    for (long long j = 0; j < Nsites; j++)
    {

        Vec3 site;

        do
        {
            site = Vec3(
                rng.Uniform(0.0, L),
                rng.Uniform(-W / 2.0, W / 2.0),
                rng.Uniform(-T / 2.0, T / 2.0));
        } while (!InsideActiveVolume(site.x, site.y, site.z, cfg));

        for (long long i = 0; i < NphotPerSite; i++)
        {
            PhotonResult res = PropagateOnePhoton(rng, site, (int)j,  cfg);
            if (res.detected)
                nDet++;
            wr.Fill(res);
        }
    }

    wr.Close();

    double frac = (Nsites > 0 && NphotPerSite > 0) ? double(nDet) / double(Nsites * NphotPerSite) : 0.0;
    std::cout << "RunManySites: detected fraction = " << frac
              << " (Nsites=" << Nsites << ", NphotPerSite=" << NphotPerSite << ")\n";
}

