import ROOT, os

# Input file
dir  = "/Users/veronicapilloni/Desktop/EXO/"
file = "nano.root"
events = ROOT.RDataFrame("Events", os.path.join(dir, file))

# PDG IDs
pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
deltaR_threshold = 0.5

# C++ helper functions
ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

using ROOT::VecOps::RVec;

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = std::fabs(phi1 - phi2);
    if(dphi > M_PI) dphi = 2*M_PI - dphi;
    double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// Flag per ogni dark pion: 1 se → 2 dark photons → 2 muoni ciascuno, 0 altrimenti
RVec<int> IsGoodDarkPion(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkpion,
    int pdgID_darkphoton
) {
    RVec<int> flags(pdgId.size(), 0);

    // conta muoni per ogni dark photon
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    // marca i dark pioni "buoni"
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;

        std::vector<int> daughters;
        for (size_t j = 0; j < pdgId.size(); ++j) {
            if (motherIdx[j] == (int)i && pdgId[j] == pdgID_darkphoton) {
                daughters.push_back(j);
            }
        }
        if (daughters.size() != 2) continue;

        if (muon_counts[daughters[0]] == 2 && muon_counts[daughters[1]] == 2) {
            flags[i] = 1;
        }
    }

    return flags;
}

// Match reconstructed 4-muon SVs to good dark pions
RVec<int> MatchFourMuonSVsFromDarkPions(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    const RVec<float>& mu3pt,
    const RVec<float>& mu3eta,
    const RVec<float>& mu3phi,
    const RVec<float>& mu4pt,
    const RVec<float>& mu4eta,
    const RVec<float>& mu4phi,
    int pdgID_darkpion,
    int pdgID_darkphoton,
    double mu_mass,
    double dr_threshold
) {
    // Step 1: muon counts per dark photon
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    // Step 2: select good dark pions
    std::vector<int> good_darkpions;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;
        std::vector<int> daughters;
        for (size_t j = 0; j < pdgId.size(); ++j) {
            if (motherIdx[j] == (int)i && pdgId[j] == pdgID_darkphoton) {
                daughters.push_back(j);
            }
        }
        if (daughters.size() != 2) continue;
        if (muon_counts[daughters[0]] == 2 && muon_counts[daughters[1]] == 2) {
            good_darkpions.push_back(i);
        }
    }

    // Step 3: build TLorentzVector for each four-muon SV
    std::vector<TLorentzVector> fourmuSVs;
    size_t n_cands = std::min({mu1pt.size(), mu2pt.size(), mu3pt.size(), mu4pt.size()});
    for (size_t i = 0; i < n_cands; ++i) {
        TLorentzVector v1, v2, v3, v4;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        v3.SetPtEtaPhiM(mu3pt[i], mu3eta[i], mu3phi[i], mu_mass);
        v4.SetPtEtaPhiM(mu4pt[i], mu4eta[i], mu4phi[i], mu_mass);
        fourmuSVs.push_back(v1 + v2 + v3 + v4);
    }

    // Step 4: match 4μ SVs to good dark pions
    RVec<int> flags(fourmuSVs.size(), 0);
    for (size_t iSV = 0; iSV < fourmuSVs.size(); ++iSV) {
        for (int ip : good_darkpions) {
            if (deltaR(fourmuSVs[iSV].Eta(), fourmuSVs[iSV].Phi(), gp_eta[ip], gp_phi[ip]) < dr_threshold) {
                flags[iSV] = 1;
                break;
            }
        }
    }
    return flags;
}

// Compute pt of 4-muon SVs
RVec<float> ComputeFourMuonSVpt(
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    const RVec<float>& mu3pt,
    const RVec<float>& mu3eta,
    const RVec<float>& mu3phi,
    const RVec<float>& mu4pt,
    const RVec<float>& mu4eta,
    const RVec<float>& mu4phi,
    double mu_mass
) {
    size_t n_cands = std::min({mu1pt.size(), mu2pt.size(), mu3pt.size(), mu4pt.size()});
    RVec<float> pts(n_cands, 0.0);
    for (size_t i = 0; i < n_cands; ++i) {
        TLorentzVector v1, v2, v3, v4;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        v3.SetPtEtaPhiM(mu3pt[i], mu3eta[i], mu3phi[i], mu_mass);
        v4.SetPtEtaPhiM(mu4pt[i], mu4eta[i], mu4phi[i], mu_mass);
        pts[i] = (v1+v2+v3+v4).Pt();
    }
    return pts;
}

RVec<int> MatchDarkPionFlags(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<float>& gen_eta,
    const RVec<float>& gen_phi,
    const RVec<float>& mu1pt,
    const RVec<float>& mu1eta,
    const RVec<float>& mu1phi,
    const RVec<float>& mu2pt,
    const RVec<float>& mu2eta,
    const RVec<float>& mu2phi,
    const RVec<float>& mu3pt,
    const RVec<float>& mu3eta,
    const RVec<float>& mu3phi,
    const RVec<float>& mu4pt,
    const RVec<float>& mu4eta,
    const RVec<float>& mu4phi,
    int pdgID_darkpion,
    int pdgID_darkphoton,
    double mu_mass,
    double dr_threshold
) {
    RVec<int> matched_flags(pdgId.size(), 0);

    // ---- Step 1: costruisci i TLorentzVector dei fourmuonSV ----
    std::vector<TLorentzVector> fourmuSVs;
    size_t n_cands = std::min({mu1pt.size(), mu2pt.size(), mu3pt.size(), mu4pt.size()});
    for (size_t i = 0; i < n_cands; ++i) {
        TLorentzVector v1, v2, v3, v4;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        v3.SetPtEtaPhiM(mu3pt[i], mu3eta[i], mu3phi[i], mu_mass);
        v4.SetPtEtaPhiM(mu4pt[i], mu4eta[i], mu4phi[i], mu_mass);
        fourmuSVs.push_back(v1 + v2 + v3 + v4);
    }

    // ---- Step 2: trova i dark pions "buoni" ----
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == 13) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }

    std::vector<int> good_darkpions;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;

        std::vector<int> daughters;
        for (size_t j = 0; j < pdgId.size(); ++j) {
            if (motherIdx[j] == (int)i && pdgId[j] == pdgID_darkphoton) {
                daughters.push_back(j);
            }
        }
        if (daughters.size() != 2) continue;

        if (muon_counts[daughters[0]] == 2 && muon_counts[daughters[1]] == 2) {
            good_darkpions.push_back(i);
        }
    }

    // ---- Step 3: match dark pions buoni con fourmuonSV ----
    for (int ip : good_darkpions) {
        for (auto& sv : fourmuSVs) {
            if (deltaR(sv.Eta(), sv.Phi(), gen_eta[ip], gen_phi[ip]) < dr_threshold) {
                matched_flags[ip] = 1;
                break;
            }
        }
    }

    return matched_flags;
}
""")

# DataFrame pipeline
def prepare_dataframe(df_raw):
    return (
        df_raw
        .Define("isGoodDarkPion",
            f"IsGoodDarkPion(GenPart_pdgId, GenPart_genPartIdxMother, {pdgID_darkpion}, {pdgID_darkphoton})")
        .Define("goodDarkPion_pt", "GenPart_pt[isGoodDarkPion == 1]")
        .Define("darkPion_flags", 
                f"MatchDarkPionFlags(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"GenPart_eta, GenPart_phi, "
                f"fourmuonSV_mu1pt, fourmuonSV_mu1eta, fourmuonSV_mu1phi, "
                f"fourmuonSV_mu2pt, fourmuonSV_mu2eta, fourmuonSV_mu2phi, "
                f"fourmuonSV_mu3pt, fourmuonSV_mu3eta, fourmuonSV_mu3phi, "
                f"fourmuonSV_mu4pt, fourmuonSV_mu4eta, fourmuonSV_mu4phi, "
                f"{pdgID_darkpion}, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        .Define("darkPion_pt_matched", "GenPart_pt[(isGoodDarkPion==1)&&(darkPion_flags == 1)]")

        
        .Define("fourmuSV_flags",
                f"MatchFourMuonSVsFromDarkPions(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"GenPart_eta, GenPart_phi, "
                f"fourmuonSV_mu1pt, fourmuonSV_mu1eta, fourmuonSV_mu1phi, "
                f"fourmuonSV_mu2pt, fourmuonSV_mu2eta, fourmuonSV_mu2phi, "
                f"fourmuonSV_mu3pt, fourmuonSV_mu3eta, fourmuonSV_mu3phi, "
                f"fourmuonSV_mu4pt, fourmuonSV_mu4eta, fourmuonSV_mu4phi, "
                f"{pdgID_darkpion}, {pdgID_darkphoton}, 0.105, {deltaR_threshold})")
        .Define("fourmuSV_pt",
                f"ComputeFourMuonSVpt(fourmuonSV_mu1pt, fourmuonSV_mu1eta, fourmuonSV_mu1phi, "
                f"fourmuonSV_mu2pt, fourmuonSV_mu2eta, fourmuonSV_mu2phi, "
                f"fourmuonSV_mu3pt, fourmuonSV_mu3eta, fourmuonSV_mu3phi, "
                f"fourmuonSV_mu4pt, fourmuonSV_mu4eta, fourmuonSV_mu4phi, 0.105)")
        .Define("fourmuSV_pt_matched", "fourmuSV_pt[fourmuSV_flags == 1]")
        .Define("nFourMuonSV_matched", "Sum(fourmuSV_flags)")
    )

df = prepare_dataframe(events)


histo = df.Histo1D(
    ("h_all", "All dark pion into 2 dark photons into 2 mu  p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "goodDarkPion_pt"
).Clone()

histo_matched = df.Histo1D(
    ("h_matched", "Matched dark photon p_{T};p_{T} [GeV];Entries", 90, 0, 60),
    "darkPion_pt_matched"
).Clone()
histo_matched.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas()
histo.Draw()
histo_matched.Draw("SAME")

print(histo.Integral())
print(histo_matched.Integral())

eff_perc = (histo_matched.Integral()/histo.Integral()) * 100
print(f"Effciency = {eff_perc:.4f} %")

c.SaveAs("DarkPions_pt_Gen_vs_Matched.png")

eff = ROOT.TEfficiency(histo_matched, histo)
c2 = ROOT.TCanvas()
eff.Draw("AP")  # "AP" = Axis + Points
c2.Update()
c2.SaveAs("EFF_DarkPions_pt_Gen_vs_Matched.png")
