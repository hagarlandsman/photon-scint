#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <sstream>

void ReadAndFillTree(const char* infile = "random_points.csv", const char* outfile = "results.root")
{
    // Variables
    double Ntot, NdetFrac, err, NabsFrac, NescFrac, NhPMTFrac;
    double cost, fDet_total, cost_per_eff_area;
    double L, W, T;
    int useWedge;
    double wedgeLen, wedgeTipW, nScint, absLen, Rwrap;
    int wrap;
    double rPMT;

    // Open output file and create a tree
    TFile* f = new TFile(outfile, "RECREATE");
    TTree* t = new TTree("t", "Photon simulation data");

    // Define branches
    t->Branch("Ntot", &Ntot);
    t->Branch("NdetFrac", &NdetFrac);
    t->Branch("err", &err);
    t->Branch("NabsFrac", &NabsFrac);
    t->Branch("NescFrac", &NescFrac);
    t->Branch("NhPMTFrac", &NhPMTFrac);
    t->Branch("cost", &cost);
    t->Branch("fDet_total", &fDet_total);
    t->Branch("cost_per_eff_area", &cost_per_eff_area);
    t->Branch("L", &L);
    t->Branch("W", &W);
    t->Branch("T", &T);
    t->Branch("useWedge", &useWedge);
    t->Branch("wedgeLen", &wedgeLen);
    t->Branch("wedgeTipW", &wedgeTipW);
    t->Branch("nScint", &nScint);
    t->Branch("absLen", &absLen);
    t->Branch("Rwrap", &Rwrap);
    t->Branch("wrap", &wrap);
    t->Branch("rPMT", &rPMT);

    // Read the file
    std::ifstream fin(infile);
    std::string line;
    int nlines = 0;
    while (std::getline(fin, line)) {
        // Skip header lines
        if (line.find("Ntot") != std::string::npos) continue;

        std::stringstream ss(line);
        char comma;

        // Read all columns
        if ((ss >> Ntot >> comma
                  >> NdetFrac >> comma
                  >> err >> comma
                  >> NabsFrac >> comma
                  >> NescFrac >> comma
                  >> NhPMTFrac >> comma
                  >> cost >> comma
                  >> fDet_total >> comma
                  >> cost_per_eff_area >> comma
                  >> L >> comma
                  >> W >> comma
                  >> T >> comma
                  >> useWedge >> comma
                  >> wedgeLen >> comma
                  >> wedgeTipW >> comma
                  >> nScint >> comma
                  >> absLen >> comma
                  >> Rwrap >> comma
                  >> wrap >> comma
                  >> rPMT )) {
            nlines++;
                    t->Fill();


        }

        // Fill the tree
    }

    // Write and close
    t->Write();
    f->Close();

    std::cout << "Finished "<<nlines<<" lines. Reading and saving to " << outfile << std::endl;
}
