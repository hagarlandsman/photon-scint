// PrintCountsFromTree.C
#include "TFile.h"
#include "TTree.h"
#include <iostream>

void PrintCountsFromTree(const char* fn="toyOptics.root", const char* tn="tPhot")
{
    TFile f(fn, "READ");
    if (f.IsZombie()) {
        std::cout << "Cannot open file: " << fn << "\n";
        return;
    }

    auto t = (TTree*)f.Get(tn);
    if (!t) {
        std::cout << "Cannot find tree: " << tn << "\n";
        return;
    }

    auto has = [&](const char* b){ return t->GetBranch(b) != nullptr; };

    // You might have one of these depending on how you wrote the tree
    // absorbed: "absorbed"
    // escaped:  "escaped"
    // hitPMT:   "inPMT" or "hitPMT"
    // detected: "detected"
    const char* bAbs = has("absorbed") ? "absorbed" : nullptr;
    const char* bEsc = has("escaped")  ? "escaped"  : nullptr;
    const char* bHit = has("inPMT")    ? "inPMT"    : (has("hitPMT") ? "hitPMT" : nullptr);
    const char* bDet = has("detected") ? "detected" : nullptr;

    if (!bAbs) std::cout << "Missing branch: absorbed\n";
    if (!bEsc) std::cout << "Missing branch: escaped\n";
    if (!bHit) std::cout << "Missing branch: inPMT or hitPMT\n";
    if (!bDet) std::cout << "Missing branch: detected\n";

    auto countOnes = [&](const char* bname)->Long64_t {
        if (!bname) return 0;
        // assumes 0/1 per entry
        TString sel = Form("%s!=0", bname);
        return t->GetEntries(sel);
    };

    Long64_t nTot = t->GetEntries();

    Long64_t nAbs = countOnes(bAbs);
    Long64_t nEsc = countOnes(bEsc);
    Long64_t nHit = countOnes(bHit);
    Long64_t nDet = countOnes(bDet);

    std::cout << "File: " << fn << " | Tree: " << tn << "\n";
    std::cout << "Total entries: " << nTot << "\n";
    std::cout << "absorbed: "  << nAbs << "\n";
    std::cout << "escaped: "   << nEsc << "\n";
    std::cout << "hitPMT: "    << nHit << "  (branch=" << (bHit ? bHit : "none") << ")\n";
    std::cout << "detected: "  << nDet << "\n";
}
