import ROOT, os

# Input file
dir  = "/Users/veronicapilloni/Desktop/EXO/"
file = "nano.root"
events = ROOT.RDataFrame("Events", os.path.join(dir, file))

# PDG IDs
pdgID_darkphoton = 999999
deltaR_threshold = 0.5


# Declare TLorentzVector and deltaR in C++
ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
#include <algorithm>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

using ROOT::VecOps::RVec;

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = std::fabs(phi1 - phi2);
    if(dphi > M_PI) dphi = 2*M_PI - dphi;
    double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// Restituisce un flag per ogni GenPart (1 = dark photon matchato, 0 altrimenti)
RVec<int> MatchPhotonFlags(
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    const RVec<int>& pdgId,
    const RVec<float>& gp_pt,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    int pdgID_darkphoton,
    double mu_mass,
    double dr_threshold
) {
    RVec<int> matched_flags(gp_pt.size(), 0);

    std::vector<double> dimu_eta;
    std::vector<double> dimu_phi;
    size_t n_pairs = std::min(mu1pt.size(), mu2pt.size());
    for(size_t i=0; i<n_pairs; ++i){
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        TLorentzVector dimu = v1 + v2;
        dimu_eta.push_back(dimu.Eta());
        dimu_phi.push_back(dimu.Phi());
    }

    for(size_t ip=0; ip<gp_pt.size(); ++ip){
        if(pdgId[ip] != pdgID_darkphoton) continue;
        for(size_t j=0; j<dimu_eta.size(); ++j){
            if(deltaR(dimu_eta[j], dimu_phi[j], gp_eta[ip], gp_phi[ip]) < dr_threshold){
                matched_flags[ip] = 1;
                break;
            }
        }
    }
    return matched_flags;
}


""")


def prepare_dataframe(df_raw):
    return(
        df_raw
        #photons
        .Define("isDarkPhoton", f"GenPart_pdgId == {pdgID_darkphoton}")
        .Define("nDarkPhotons",  f"Sum(GenPart_pdgId == {pdgID_darkphoton})")
        .Define("darkPhoton_pt", "GenPart_pt[isDarkPhoton]")
        .Define("darkPhoton_eta", "GenPart_eta[isDarkPhoton]")
        .Define("darkPhoton_phi", "GenPart_phi[isDarkPhoton]")
        df = df.Define("darkPhoton_flags", f"MatchPhotonFlags(muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        df = df.Define("darkPhoton_pt_matched", "GenPart_pt[isDarkPhoton&&(darkPhoton_flags ==1)]")

)
        
df = prepare_dataframe(events)

histo = df.Histo1D(
    ("h_all", "All dark photon p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "darkPhoton_pt"
).Clone()

histo_matched = df.Histo1D(
    ("h_matched", "Matched dark photon p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "darkPhoton_pt_matched"
).Clone()
histo_matched.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas()
histo.Draw()
histo_matched.Draw("SAME")

print(histo.Integral())
print(histo_matched.Integral())
c.SaveAs("DarkPhotons_pt_Gen_vs_Matched.png")

eff = ROOT.TEfficiency(histo_matched, histo)
c2 = ROOT.TCanvas()
eff.Draw("AP")  # "AP" = Axis + Points
c2.Update()
c2.SaveAs("EFF_DarkPhotons_pt_Gen_vs_Matched.png")
