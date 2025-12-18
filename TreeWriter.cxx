#include "TFile.h"
#include "TTree.h"
#include "OpticsConfig.h"
#include "TreeWriter.h"
#include <iostream>


TreeWriter::TreeWriter(const char *outFile, const OpticsConfig &cfg)
        : cfg_(cfg), savePath_(cfg.savePath)

    {
        fout_ = new TFile(outFile, "RECREATE");
        t_ = new TTree("tPhot", "Toy optical photon transport");

        t_->Branch("site_number", &r_.site_number, "site_number/I");

        t_->Branch("x0", &r_.x0, "x0/D");
        t_->Branch("y0", &r_.y0, "y0/D");
        t_->Branch("z0", &r_.z0, "z0/D");

        t_->Branch("xf", &r_.xf, "xf/D");
        t_->Branch("yf", &r_.yf, "yf/D");
        t_->Branch("zf", &r_.zf, "zf/D");

        t_->Branch("epsRaysPMT", &r_.epsRaysPMT, "epsRaysPMT/D");

        t_->Branch("endPlane", &r_.endPlane, "endPlane/I");
        t_->Branch("pmt_side", &r_.pmt_side, "pmt_side/I");
        t_->Branch("inPMT", &r_.inPMT, "inPMT/I");
        t_->Branch("detected", &r_.detected, "detected/I");
        t_->Branch("absorbed", &r_.absorbed, "absorbed/I");
        t_->Branch("escaped", &r_.escaped, "escaped/I");
        t_->Branch("reachedEnd", &r_.reachedEnd, "reachedEnd/I");

        t_->Branch("nBounces", &r_.nBounces, "nBounces/I");
        t_->Branch("path", &r_.path, "path/D");

        if (savePath_)
        {
            t_->Branch("xPath", &r_.xPath);
            t_->Branch("yPath", &r_.yPath);
            t_->Branch("zPath", &r_.zPath);
        }
                // 2) Add config branches (repeated per entry)
        // geometry + wrap live in cfg now (if you do your refactor)
        t_->Branch("cfg_L", &cfg_L_, "cfg_L/D");
        t_->Branch("cfg_W", &cfg_W_, "cfg_W/D");
        t_->Branch("cfg_T", &cfg_T_, "cfg_T/D");
        t_->Branch("cfg_wrap", &cfg_wrap_, "cfg_wrap/I");

        t_->Branch("cfg_nScint", &cfg_nScint_, "cfg_nScint/D");
        t_->Branch("cfg_nOut", &cfg_nOut_, "cfg_nOut/D");
        t_->Branch("cfg_absLen", &cfg_absLen_, "cfg_absLen/D");
        t_->Branch("cfg_Rwrap", &cfg_Rwrap_, "cfg_Rwrap/D");
        t_->Branch("cfg_maxSteps", &cfg_maxSteps_, "cfg_maxSteps/I");

        t_->Branch("cfg_rPMT", &cfg_rPMT_, "cfg_rPMT/D");
        t_->Branch("cfg_epsCouple", &cfg_epsCouple_, "cfg_epsCouple/D");
        t_->Branch("cfg_pde", &cfg_pde_, "cfg_pde/D");

        t_->Branch("cfg_eps0", &cfg_eps0_, "cfg_eps0/D");
        t_->Branch("cfg_lambdaC", &cfg_lambdaC_, "cfg_lambdaC/D");

        t_->Branch("cfg_savePath", &cfg_savePath_, "cfg_savePath/I");

        t_->Branch("cfg_useWedge", &cfg_useWedge_, "cfg_useWedge/I");
        t_->Branch("cfg_wedgeLen", &cfg_wedgeLen_, "cfg_wedgeLen/D");
        t_->Branch("cfg_wedgeTipW", &cfg_wedgeTipW_, "cfg_wedgeTipW/D");

        // freeze the config values once
        CopyCfgToBranchScalars_();


    }

    void TreeWriter::Fill(const PhotonResult &in)
    {
        r_ = in;
              // ensure cfg scalars are set (cheap; safe if you ever mutate cfg_)
        // you can remove this call if you guarantee cfg_ never changes
        CopyCfgToBranchScalars_();
        if (!savePath_)
        {
            r_.xPath.clear();
            r_.yPath.clear();
            r_.zPath.clear();
        }
        t_->Fill();
    }

    void TreeWriter::Close()
    {
        fout_->cd();
        t_->Write();
        fout_->Close();
        delete fout_;
        fout_ = nullptr;
        t_ = nullptr;
    }

    void TreeWriter::CopyCfgToBranchScalars_()
    {
        cfg_L_ = cfg_.L;
        cfg_W_ = cfg_.W;
        cfg_T_ = cfg_.T;
        cfg_wrap_ = cfg_.wrap;

        cfg_nScint_ = cfg_.nScint;
        cfg_nOut_ = cfg_.nOut;
        cfg_absLen_ = cfg_.absLen;
        cfg_Rwrap_ = cfg_.Rwrap;
        cfg_maxSteps_ = cfg_.maxSteps;

        cfg_rPMT_ = cfg_.rPMT;
        cfg_epsCouple_ = cfg_.epsCouple;
        cfg_pde_ = cfg_.pde;

        cfg_eps0_ = cfg_.eps0;
        cfg_lambdaC_ = cfg_.lambdaC;

        cfg_savePath_ = cfg_.savePath ? 1 : 0;

        cfg_useWedge_ = cfg_.useWedge ? 1 : 0;
        cfg_wedgeLen_ = cfg_.wedgeLen;
        cfg_wedgeTipW_ = cfg_.wedgeTipW;
    }


