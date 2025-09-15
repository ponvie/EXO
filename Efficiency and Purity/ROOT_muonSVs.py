import ROOT, os

# Input file
dir  = "/Users/veronicapilloni/Desktop/EXO/"
#file = "nano.root"
file = "data_0.root"
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

RVec<int> CountMuonPhotons(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkphoton
) {
    RVec<int> muon_counts(pdgId.size(), 0);

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    return muon_counts;
}



RVec<int> MatchMuonSVsFrom2MuDarkPhotons(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<float>& gp_pt,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    int pdgID_darkphoton,
    double mu_mass,
    double dr_threshold
) {
    // Step 1: conta muoni figli dei dark photon
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    // Step 2: costruisci i 4-vettori dei muonSV
    std::vector<TLorentzVector> dimuons;
    size_t n_pairs = std::min(mu1pt.size(), mu2pt.size());
    for (size_t i = 0; i < n_pairs; ++i) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        dimuons.push_back(v1 + v2);
    }

    // Step 3: crea flag per ogni muonSV
    RVec<int> flags(dimuons.size(), 0);
    for (size_t iSV = 0; iSV < dimuons.size(); ++iSV) {
        for (size_t ip = 0; ip < gp_pt.size(); ++ip) {
            if (pdgId[ip] != pdgID_darkphoton) continue;
            if (muon_counts[ip] != 2) continue; // solo dark photon → 2 muoni

            if (deltaR(dimuons[iSV].Eta(), dimuons[iSV].Phi(), gp_eta[ip], gp_phi[ip]) < dr_threshold) {
                flags[iSV] = 1;
                break; // basta un match
            }
        }
    }

    return flags; // mask da usare su muonSV
}

RVec<float> ComputeDiMuonPt(
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    double mu_mass
) {
    size_t n_pairs = std::min(mu1pt.size(), mu2pt.size());
    RVec<float> dimuon_pt(n_pairs, 0.0);

    for (size_t i = 0; i < n_pairs; ++i) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        TLorentzVector dimu = v1 + v2;
        dimuon_pt[i] = dimu.Pt();
    }

    return dimuon_pt;
}

int CountMatchedMuonSVsFrom2MuDarkPhotons(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<float>& gp_pt,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    int pdgID_darkphoton,
    double mu_mass,
    double dr_threshold
) {
    // Step 1: conta muoni figli dei dark photon
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    // Step 2: costruisci i 4-vettori dei muonSV
    std::vector<TLorentzVector> dimuons;
    size_t n_pairs = std::min(mu1pt.size(), mu2pt.size());
    for (size_t i = 0; i < n_pairs; ++i) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        dimuons.push_back(v1 + v2);
    }

    // Step 3: conta quanti muonSV sono matchati a dark photon che decadono in 2 muoni
    int matched_count = 0;
    for (size_t ip = 0; ip < gp_pt.size(); ++ip) {
        if (pdgId[ip] != pdgID_darkphoton) continue;
        if (muon_counts[ip] != 2) continue; // solo dark photon → 2 muoni

        for (auto& dimu : dimuons) {
            if (deltaR(dimu.Eta(), dimu.Phi(), gp_eta[ip], gp_phi[ip]) < dr_threshold) {
                matched_count++;
                break; // basta un match per dark photon
            }
        }
    }

    return matched_count;
}


""")


def prepare_dataframe(df_raw):
    return(
        df_raw
        #photons
        .Define("muon_counts", f"CountMuonPhotons(GenPart_pdgId, GenPart_genPartIdxMother, {pdgID_darkphoton})")
        
        #.Define("isDarkPhoton", f"GenPart_pdgId == {pdgID_darkphoton}")
        .Define("isDarkPhoton_into2mu", f"(GenPart_pdgId == {pdgID_darkphoton}) && (muon_counts == 2) ")
        
        #.Define("nDarkPhotons",  f"Sum(GenPart_pdgId == {pdgID_darkphoton})")
        .Define("nDarkPhoton_into2mu", f"Sum(GenPart_pdgId == {pdgID_darkphoton}) && (muon_counts == 2)")
        
        #.Define("darkPhoton_pt", "GenPart_pt[isDarkPhoton]")
        #.Define("darkPhoton_eta", "GenPart_eta[isDarkPhoton]")
        #.Define("darkPhoton_phi", "GenPart_phi[isDarkPhoton]")
        
        .Define("darkPhoton_pt", "GenPart_pt[isDarkPhoton_into2mu]")
        .Define("darkPhoton_eta", "GenPart_eta[isDarkPhoton_into2mu]")
        .Define("darkPhoton_phi", "GenPart_phi[isDarkPhoton_into2mu]")
        
        .Define("darkPhoton_flags", f"MatchPhotonFlags(muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        #.Define("darkPhoton_pt_matched", "GenPart_pt[isDarkPhoton&&(darkPhoton_flags ==1)]")
        .Define("darkPhoton_pt_matched", "GenPart_pt[isDarkPhoton_into2mu&& (darkPhoton_flags == 1)]")
        
        
        .Define("nMatchedMuonSVs_from2muDarkPhoton", f"CountMatchedMuonSVsFrom2MuDarkPhotons(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        .Define("muonSVs_flags", f"MatchMuonSVsFrom2MuDarkPhotons(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        .Define("muonSVs_pt", f"ComputeDiMuonPt(muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, 0.105)")        
        .Define("muonSVs_pt_matched", "muonSVs_pt[muonSVs_flags==1]")
        .Define("nmuonSV_matched",  "Sum(muonSVs_flags==1)")
    
)         
df = prepare_dataframe(events)

#histo = df.Histo1D(tuple, variable)
#tuple = ("random_name", "histogram_title", number_of_bins, xmin, xmax)

histo = df.Histo1D(
    ("h_all", "All dark photon p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "darkPhoton_pt"
).Clone()

#Efficiency
histo_matched = df.Histo1D(
    ("h_matched", "Matched dark photon p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "darkPhoton_pt_matched"
).Clone()
histo_matched.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas()
histo.Draw()
histo_matched.Draw("SAME")
c.SaveAs("DarkPhotons_pt_Gen_vs_Matched.png")
print(histo.Integral())
print(histo_matched.Integral())

eff_perc = (histo_matched.Integral()/histo.Integral()) * 100
print(f"Effciency = {eff_perc:.4f} %")
eff = ROOT.TEfficiency(histo_matched, histo)
c2 = ROOT.TCanvas()
eff.Draw("AP")  # "AP" = Axis + Points
c2.Update()
c2.SaveAs("EFF_DarkPhotons_pt_Gen_vs_Matched.png")


nSV_histo = df.Histo1D(
    ("h_all", "All muonSVs;Entries", 10, 0.50, 10.5),
    "nmuonSV"
).Clone()
nSV_matched_histo = df.Histo1D(
    ("h_matched", "Matched muonSVs;Entries", 10, 0.50, 10.5),
    "nmuonSV_matched"
).Clone()

c3 = ROOT.TCanvas()
nSV_histo.Draw()
nSV_matched_histo.Draw("SAME")
nSV_matched_histo.SetLineColor(ROOT.kRed)
nSV_matched_histo.SetLineStyle(2)

c3.Update()
c3.SaveAs("nMuonSVs_distr.png")

print("nMuonSVs:", nSV_histo.Integral(0, 91))
print("nMuonSVs matched:", nSV_matched_histo.Integral(0, 91))

#Purity
print("total muonSVs:", df.Sum("nmuonSV").GetValue())
print("muonsSVs matched to dark photons that decay into two muons:", df.Sum("nMatchedMuonSVs_from2muDarkPhoton").GetValue())

histo_muonSVs = df.Histo1D(
    ("h_all", "All muonSVs p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "muonSVs_pt"
).Clone()

histo_muonSVs_matched = df.Histo1D(
    ("h_matched", "Matched muonSVs p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "muonSVs_pt_matched"
).Clone()
histo_muonSVs_matched.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas()
histo_muonSVs.Draw()
histo_muonSVs_matched.Draw("SAME")
c.SaveAs("MuonSVs_pt.png")

print(histo_muonSVs.Integral())
print(histo_muonSVs_matched.Integral())

purity = (histo_muonSVs_matched.Integral()/histo_muonSVs.Integral()) * 100
print(f"Purity = {purity:.4f} %")

pur=ROOT.TEfficiency(histo_muonSVs_matched,histo_muonSVs)
c3 = ROOT.TCanvas()
pur.Draw("AP")  # "AP" = Axis + Points
c3.Update()
c3.SaveAs("PUR_muonSVs_pt_vs_Matched.png")
