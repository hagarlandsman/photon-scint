// TreeWriter.h
#pragma once

#include "TFile.h"
#include "TTree.h"

#include "OpticsConfig.h"   // defines PhotonResult (incl. xPath/yPath/zPath)

class TreeWriter
{
public:
    TreeWriter(const char *outFile, const OpticsConfig &cfg);

 //   ~TreeWriter();

    void Fill(const PhotonResult &in);
    void Close();

private:
    TFile *fout_ = nullptr;
    TTree *t_ = nullptr;
    PhotonResult r_;
    bool savePath_ = false;
    bool closed_ = false;
    void CopyCfgToBranchScalars_();

    OpticsConfig cfg_;

    // scalars used as branch addresses
    double cfg_L_ = 0, cfg_W_ = 0, cfg_T_ = 0;
    int cfg_wrap_ = 1;

    double cfg_nScint_ = 0, cfg_nOut_ = 0;
    double cfg_absLen_ = 0, cfg_Rwrap_ = 0;
    int cfg_maxSteps_ = 0;

    double cfg_rPMT_ = 0, cfg_epsCouple_ = 0, cfg_pde_ = 0;

    double cfg_eps0_ = 0, cfg_lambdaC_ = 0;
    int cfg_savePath_ = 0;

    int cfg_useWedge_ = 0;
    double cfg_wedgeLen_ = 0, cfg_wedgeTipW_ = 0;
};
