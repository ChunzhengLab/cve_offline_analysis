// lambda_correlation_ratio.C
// -----------------------------------------------------------------------------
// Compute the ratio
//      R(Δφ) = ⟨cos(φ_Λ1 − φ_p2)⟩ / ⟨cos(φ_Λ1 − φ_Λ2)⟩
// as a function of Δφ(Λ1,Λ2) = φ_Λ1 − φ_Λ2, **requiring pT(Λ1) > 1 GeV/c and
// pT(Λ2) > 1 GeV/c**.
//
//  - Reads two TH1D spectra (Λ1 and Λ2 pT) from an input ROOT file.
//  - Generates isotropic Λ2 → p + π− decays; the proton is boosted to the lab.
//  - Uses TProfile to accumulate the numerator and denominator separately,
//    then forms their bin‑by‑bin ratio.
//  - Stores pNum, pDen, and hRatio in an output ROOT file, and saves a PNG.
//
// Compile (ACLiC):
//      root -l 'lambda_correlation_ratio.C++'
// Call (interpreted):
//      root -l -q 'lambda_correlation_ratio.C("in.root","hpt1","hpt2",1e6)'
// -----------------------------------------------------------------------------

#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

// ----- PDG masses (GeV/c²) ---------------------------------------------------
constexpr double M_LAMBDA = 1.115683;
constexpr double M_PROTON = 0.938272;
constexpr double M_PION   = 0.139570;

// ----------------------------------------------------------------------------
// Helper: isotropic Λ → p + π in Λ rest frame, then boost to lab
// ----------------------------------------------------------------------------
inline TLorentzVector ProtonFromLambda(const TLorentzVector& lambda, TRandom3& rng)
{
    const double p_mag = std::sqrt((M_LAMBDA*M_LAMBDA - std::pow(M_PROTON+M_PION,2)) *
                                   (M_LAMBDA*M_LAMBDA - std::pow(M_PROTON-M_PION,2)) ) /
                         (2.0 * M_LAMBDA);

    const double costh = 2.0 * rng.Rndm() - 1.0;
    const double sinth = std::sqrt(1.0 - costh*costh);
    const double phiRF = rng.Uniform(0.0, 2.0*TMath::Pi());

    TVector3 pRF(p_mag * sinth * std::cos(phiRF),
                 p_mag * sinth * std::sin(phiRF),
                 p_mag * costh);

    TLorentzVector p4(pRF, std::sqrt(pRF.Mag2() + M_PROTON*M_PROTON));
    p4.Boost(lambda.BoostVector());
    return p4;
}

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void lambda_correlation_ratio(const char* inFile  = "/Users/wangchunzheng/works/Experiments/cve/src_fraction_fit/src/pt_dist.root",
                              const char* h1Name = "h_pt_lambda_18q_cent3_eff_col",
                              const char* h2Name = "h_pt_lambda_18q_cent3_eff_col",
                              Long64_t    nEv    = 10000000,
                              const char* outFile= "lambda_corr_ratio.root")
{
    // -- open input ---------------------------------------------------------
    TFile fin(inFile, "READ");
    if (fin.IsZombie()) {
        std::cerr << "[ERROR] Cannot open " << inFile << "\n";
        return;
    }

    TH1D* hpt1 = nullptr;
    TH1D* hpt2 = nullptr;
    fin.GetObject(h1Name, hpt1);
    fin.GetObject(h2Name, hpt2);
    if (!hpt1 || !hpt2) {
        std::cerr << "[ERROR] Missing histogram(s) '" << h1Name << "' or '" << h2Name << "'\n";
        return;
    }

    // -- book profiles ------------------------------------------------------
    const int nBins = 64;
    TProfile pNum("pNum", "#LTcos(#phi_{#Lambda1}-#phi_{p2})#GT;#Delta#phi(#Lambda_{1},#Lambda_{2});#LTcos#GT", nBins,0.,TMath::Pi());
    TProfile pDen("pDen", "#LTcos(#phi_{#Lambda1}-#phi_{#Lambda2})#GT;#Delta#phi(#Lambda_{1},#Lambda_{2});#LTcos#GT", nBins,0.,TMath::Pi());

    // -- random generator ---------------------------------------------------
    TRandom3 rng(0);

    // -- event loop ---------------------------------------------------------
    for (Long64_t ie = 0; ie < nEv; ++ie) {
        const double pt1  = hpt1->GetRandom();
        const double pt2  = hpt2->GetRandom();

        // ---- pT cut -------------------------------------------------------
        if (pt1 <= 1.0 || pt2 <= 1.0) continue; // both Λ must satisfy pT > 1 GeV

        const double phi1 = rng.Uniform(0.0, 2.0*TMath::Pi());
        const double phi2 = rng.Uniform(0.0, 2.0*TMath::Pi());

        TLorentzVector lam1, lam2;
        lam1.SetPtEtaPhiM(pt1, 0.0, phi1, M_LAMBDA);
        lam2.SetPtEtaPhiM(pt2, 0.0, phi2, M_LAMBDA);

        // decay Λ₂ → p + π and get proton φ in lab
        const double phi_p2 = ProtonFromLambda(lam2, rng).Phi();

        const double dPhi = TVector2::Phi_mpi_pi(phi1 - phi2);

        pNum.Fill(dPhi, std::cos(TVector2::Phi_mpi_pi(phi1 - phi_p2)));
        pDen.Fill(dPhi, std::cos(dPhi));
    }

    // -- build ratio histogram ---------------------------------------------
    TH1D hRatio("hRatio","R;#Delta#phi(#Lambda_{1},#Lambda_{2});R",nBins,0,TMath::Pi());

    for (int ib = 1; ib <= nBins; ++ib) {
        const double num = pNum.GetBinContent(ib);
        const double den = pDen.GetBinContent(ib);
        if (den != 0.0) hRatio.SetBinContent(ib, num / den);
    }

    // -- write outputs ------------------------------------------------------
    TFile fout(outFile, "RECREATE");
    pNum.Write();
    pDen.Write();
    hRatio.Write();
    fout.Close();

    // -- quick look ---------------------------------------------------------
    TCanvas c("c", "Ratio vs Δφ (pT > 1 GeV/c)", 800, 500);
    hRatio.SetMarkerStyle(20);
    hRatio.SetMarkerSize(0.8);
    hRatio.Draw("PE1");
    c.SaveAs("lambda_corr_ratio_ptgt1.png");

    std::cout << "[INFO] Finished with pT > 1 GeV/c cut. Results saved to '" << outFile << "'.\n";
}
