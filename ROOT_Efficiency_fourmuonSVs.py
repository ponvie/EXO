import ROOT, os

# --- Input file ------------------------------------------------------------
dir  = "/Users/veronicapilloni/Desktop/EXO/"
file = "nano.root"
events = ROOT.RDataFrame("Events", os.path.join(dir, file))

# --- PDG IDs and parameters ------------------------------------------------
pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
deltaR_threshold = 0.5
muon_mass        = 0.105

# --- C++ helper functions --------------------------------------------------
ROOT.gInterpreter.Declare("""
// Needed headers
#include <vector>
#include <cmath>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
using ROOT::VecOps::RVec;

// --------------------------------------------------------------------------
// Utility: deltaR between two directions
double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2*M_PI - dphi;
    double deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// --------------------------------------------------------------------------
// Utility: get photon daughters of a given dark pion
std::vector<int> GetPhotonDaughters(
    size_t pion_idx,
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkphoton
) {
    std::vector<int> daughters;
    for (size_t j = 0; j < pdgId.size(); ++j) {
        if (motherIdx[j] == (int)pion_idx && pdgId[j] == pdgID_darkphoton) {
            daughters.push_back(j);
        }
    }
    return daughters;
}

// --------------------------------------------------------------------------
// Count how many dark photons each dark pion decays into
RVec<int> CountPhotonPions(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkpion,
    int pdgID_darkphoton
) {
    RVec<int> photon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] == pdgID_darkphoton) {
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkpion) {
                photon_counts[mom] += 1;
            }
        }
    }
    return photon_counts;
}

// --------------------------------------------------------------------------
// Count how many muons each dark photon decays into
RVec<int> CountMuonPhotons(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkphoton
) {
    RVec<int> muon_counts(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == 13) {  // muon
            int mom = motherIdx[i];
            if (mom >= 0 && mom < (int)pdgId.size() && pdgId[mom] == pdgID_darkphoton) {
                muon_counts[mom] += 1;
            }
        }
    }
    return muon_counts;
}

// --------------------------------------------------------------------------
// Match dark photon to dimuon SV
RVec<int> MatchPhotonToMuonSVFlag(
    const RVec<int>& pdgId,
    const RVec<int>& muon_counts,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    const RVec<float>& mu1pt, const RVec<float>& mu1eta, const RVec<float>& mu1phi,
    const RVec<float>& mu2pt, const RVec<float>& mu2eta, const RVec<float>& mu2phi,
    int pdgID_darkphoton, double mu_mass, double dr_threshold
) {
    RVec<int> matched_flags(pdgId.size(), 0);

    size_t n_vertices = std::min(mu1pt.size(), mu2pt.size());
    std::vector<TLorentzVector> dimuons;
    for (size_t i = 0; i < n_vertices; ++i) {
        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        dimuons.push_back(v1+v2);
    }

    for (size_t ip = 0; ip < pdgId.size(); ++ip) {
        if (pdgId[ip] != pdgID_darkphoton) continue;
        if (muon_counts[ip] != 2) continue; // require photon → 2 muons

        for (auto& vtx : dimuons) {
            if (deltaR(vtx.Eta(), vtx.Phi(), gp_eta[ip], gp_phi[ip]) < dr_threshold) {
                matched_flags[ip] = 1;
                break;
            }
        }
    }
    return matched_flags;
}

// --------------------------------------------------------------------------
// Flag: dark pion -> 2 photons, each photon -> 2 muons
RVec<int> IsDarkPion_in2Photons_into2Muons(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<int>& photon_counts,
    const RVec<int>& muon_counts,
    int pdgID_darkpion,
    int pdgID_darkphoton
) {
    RVec<int> flags(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;
        auto daughters = GetPhotonDaughters(i, pdgId, motherIdx, pdgID_darkphoton);
        if (daughters.size() == 2 &&
            muon_counts[daughters[0]] == 2 &&
            muon_counts[daughters[1]] == 2) {
            flags[i] = 1;
        }
    }
    return flags;
}

// --------------------------------------------------------------------------
// Flag: pion with 2 photons, requiring N of them matched to muonSVs
RVec<int> IsDarkPionPhotonMuonSVMatch(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<int>& photon_matched_flags,
    int pdgID_darkpion,
    int pdgID_darkphoton,
    int required_matches
) {
    RVec<int> flags(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;
        auto daughters = GetPhotonDaughters(i, pdgId, motherIdx, pdgID_darkphoton);
        if (daughters.size() != 2) continue;
        int nMatched = photon_matched_flags[daughters[0]] + photon_matched_flags[daughters[1]];
        if ((required_matches == 1 && nMatched >= 1) ||
            (required_matches == 2 && nMatched == 2)) {
            flags[i] = 1;
        }
    }
    return flags;
}

// --------------------------------------------------------------------------
// Match dark pions to four-muon secondary vertices (4muSV)
RVec<int> MatchPionFlags(
    const RVec<float>& mu1pt, const RVec<float>& mu1eta, const RVec<float>& mu1phi,
    const RVec<float>& mu2pt, const RVec<float>& mu2eta, const RVec<float>& mu2phi,
    const RVec<float>& mu3pt, const RVec<float>& mu3eta, const RVec<float>& mu3phi,
    const RVec<float>& mu4pt, const RVec<float>& mu4eta, const RVec<float>& mu4phi,
    const RVec<int>& pdgId,
    const RVec<float>& gp_pt, const RVec<float>& gp_eta, const RVec<float>& gp_phi,
    int pdgID_darkpion, double mu_mass, double dr_threshold
) {
    RVec<int> matched_flags(gp_pt.size(), 0);
    std::vector<TLorentzVector> fourmuons;
    size_t n_vertices = std::min({mu1pt.size(), mu2pt.size(), mu3pt.size(), mu4pt.size()});
    for (size_t i = 0; i < n_vertices; ++i) {
        TLorentzVector v1, v2, v3, v4;
        v1.SetPtEtaPhiM(mu1pt[i], mu1eta[i], mu1phi[i], mu_mass);
        v2.SetPtEtaPhiM(mu2pt[i], mu2eta[i], mu2phi[i], mu_mass);
        v3.SetPtEtaPhiM(mu3pt[i], mu3eta[i], mu3phi[i], mu_mass);
        v4.SetPtEtaPhiM(mu4pt[i], mu4eta[i], mu4phi[i], mu_mass);
        fourmuons.push_back(v1+v2+v3+v4);
    }
    for (size_t ip = 0; ip < gp_pt.size(); ++ip) {
        if (pdgId[ip] != pdgID_darkpion) continue;
        for (auto& vtx : fourmuons) {
            if (deltaR(vtx.Eta(), vtx.Phi(), gp_eta[ip], gp_phi[ip]) < dr_threshold) {
                matched_flags[ip] = 1;
                break;
            }
        }
    }
    return matched_flags;
}

// --------------------------------------------------------------------------
// Combined: pion matched to 4muonSV AND both photons matched to muonSVs
RVec<int> IsDarkPionFourMuonSVAndTwoPhotonMuonSVs(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<int>& photon_matched_flags,
    const RVec<int>& darkPion_flags,
    int pdgID_darkpion,
    int pdgID_darkphoton
) {
    RVec<int> flags(pdgId.size(), 0);
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] != pdgID_darkpion) continue;
        if (darkPion_flags[i] != 1) continue; // require 4muonSV match
        auto daughters = GetPhotonDaughters(i, pdgId, motherIdx, pdgID_darkphoton);
        if (daughters.size() == 2 &&
            photon_matched_flags[daughters[0]] == 1 &&
            photon_matched_flags[daughters[1]] == 1) {
            flags[i] = 1;
        }
    }
    return flags;
}
""")

# --- Define analysis pipeline with RDataFrame ------------------------------
def prepare_dataframe(df_raw):
    return (
        df_raw
        # Gen-level classification
        .Define("photon_counts",
                f"CountPhotonPions(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"{pdgID_darkpion}, {pdgID_darkphoton})")
        .Define("muon_counts",
                f"CountMuonPhotons(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"{pdgID_darkphoton})")

        # Flags
        .Define("isDarkPion_in2darkPhotons",
                f"(GenPart_pdgId == {pdgID_darkpion}) && (photon_counts == 2)")
        .Define("isDarkPhoton_into2muons",
                f"(GenPart_pdgId == {pdgID_darkphoton}) && (muon_counts == 2)")

        # Photon → muonSV matching
        .Define("photon_matched_flags",
               f"MatchPhotonToMuonSVFlag(GenPart_pdgId, muon_counts, "
               f"GenPart_eta, GenPart_phi, "
               f"muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi, "
               f"muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi, "
               f"{pdgID_darkphoton}, {muon_mass}, {deltaR_threshold})")
        #----------------------------------------------------------------------------------#        
        # DENOMINATOR: pion → 2 photons → 2 muons
        .Define("isDarkPion_in2darkPhotons_into2muons",
                f"IsDarkPion_in2Photons_into2Muons(GenPart_pdgId, "
                f"GenPart_genPartIdxMother, photon_counts, muon_counts, "
                f"{pdgID_darkpion}, {pdgID_darkphoton})")
        .Define("darkPion_DEN_pt", "GenPart_pt[isDarkPion_in2darkPhotons_into2muons==1]")
        #----------------------------------------------------------------------------------#        
        # NUMERATOR 1: ≥1 photon matched
        .Define("darkPion_in2darkPhotons_1muonSV",
                f"IsDarkPionPhotonMuonSVMatch(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"photon_matched_flags, {pdgID_darkpion}, {pdgID_darkphoton}, 1)")
        #.Define("darkPion_NUM1_pt", "GenPart_pt[darkPion_in2darkPhotons_1muonSV==1]")
        .Define("darkPion_NUM1_pt",
        "GenPart_pt[(isDarkPion_in2darkPhotons_into2muons==1) && (darkPion_in2darkPhotons_1muonSV==1)]")
        #----------------------------------------------------------------------------------#        
        # NUMERATOR 2: 2 photons matched
        .Define("darkPion_in2darkPhotons_2muonSV",
                f"IsDarkPionPhotonMuonSVMatch(GenPart_pdgId, GenPart_genPartIdxMother, "
                f"photon_matched_flags, {pdgID_darkpion}, {pdgID_darkphoton}, 2)")
        #.Define("darkPion_NUM2_pt", "GenPart_pt[darkPion_in2darkPhotons_2muonSV==1]")
        .Define("darkPion_NUM2_pt",
        "GenPart_pt[(isDarkPion_in2darkPhotons_into2muons==1) && (darkPion_in2darkPhotons_2muonSV==1)]")
        #----------------------------------------------------------------------------------#        
        # NUMERATOR 3: pion matched to 4muonSV
        .Define("darkPion_into4muonSV_flags",
                f"MatchPionFlags(fourmuonSV_mu1pt,fourmuonSV_mu1eta,fourmuonSV_mu1phi,"
                f"fourmuonSV_mu2pt,fourmuonSV_mu2eta,fourmuonSV_mu2phi,"
                f"fourmuonSV_mu3pt,fourmuonSV_mu3eta,fourmuonSV_mu3phi,"
                f"fourmuonSV_mu4pt,fourmuonSV_mu4eta,fourmuonSV_mu4phi,"
                f"GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi,"
                f"{pdgID_darkpion}, {muon_mass}, {deltaR_threshold})")
        #.Define("darkPion_NUM3_pt",
                #"GenPart_pt[(isDarkPion_in2darkPhotons_into2muons==1)&&(darkPion_into4muonSV_flags == 1)]")
        .Define("darkPion_NUM3_pt",
                "GenPart_pt[(isDarkPion_in2darkPhotons_into2muons==1) && (darkPion_into4muonSV_flags==1)]")
        #----------------------------------------------------------------------------------#        
        # NUMERATOR 4: pion matched to 4muonSV + 2 photons matched
        .Define("isDarkPion_FourMuonSV_and_TwoPhotonMuonSVs",
                f"IsDarkPionFourMuonSVAndTwoPhotonMuonSVs("
                f"GenPart_pdgId, GenPart_genPartIdxMother, "
                f"photon_matched_flags, darkPion_into4muonSV_flags, "
                f"{pdgID_darkpion}, {pdgID_darkphoton})")
        #.Define("darkPion_NUM4_pt", "GenPart_pt[isDarkPion_FourMuonSV_and_TwoPhotonMuonSVs==1]")
        .Define("darkPion_NUM4_pt",
                "GenPart_pt[(isDarkPion_in2darkPhotons_into2muons==1) && (isDarkPion_FourMuonSV_and_TwoPhotonMuonSVs==1)]")
        #----------------------------------------------------------------------------------#        

        # Save eta, phi
        .Define("darkPion_eta", "GenPart_eta[isDarkPion_in2darkPhotons_into2muons==1]")
        .Define("darkPion_phi", "GenPart_phi[isDarkPion_in2darkPhotons_into2muons==1]")
    )

# --- Prepare the dataframe -------------------------------------------------
df = prepare_dataframe(events)


# Counters
df = (df
      .Define("nPions_2ph",       "Sum(isDarkPion_in2darkPhotons)")
      .Define("nPions_2ph_2mu", "Sum(isDarkPion_in2darkPhotons_into2muons)")
      .Define("nPions_2ph_1muSV", "Sum(darkPion_in2darkPhotons_1muonSV)")
      .Define("nPions_2ph_2muSV",   "Sum(darkPion_in2darkPhotons_2muonSV)")
      .Define("nPions_4muSV", "Sum(darkPion_into4muonSV_flags)")
      .Define("nPions_4muSV_2ph_2muSV", "Sum(isDarkPion_FourMuonSV_and_TwoPhotonMuonSVs)") 
     )
# --- Global counts ---------------------------------------------------------
total_pions_2ph      = df.Sum("nPions_2ph").GetValue()
total_pions_2ph2mu   = df.Sum("nPions_2ph_2mu").GetValue()
total_pions_2ph_1muSV = df.Sum("nPions_2ph_1muSV").GetValue()
total_pions_2ph_2muSV = df.Sum("nPions_2ph_2muSV").GetValue()
total_pions_4muonSV = df.Sum("nPions_4muSV").GetValue()
total_pions_4muSV_2ph_2muSV = df.Sum("nPions_4muSV_2ph_2muSV").GetValue()


print("Total dark pions with 2 dark photons:", total_pions_2ph)
print("Total dark pions with 2 dark photons each one decaying into 2 muons:", total_pions_2ph2mu )
print("Total dark pions with 2 dark photons, 1 matching to a muonSVs:", total_pions_2ph_1muSV)
print("Total dark pions with 2 dark photons, 2 (both) matching to a muonSVs:", total_pions_2ph_2muSV)
print("Total dark pions matched to a 4-muonSV:", total_pions_4muonSV )
print("Total dark pions matched to a 4-muonSV, with 2 dark photons matching 2 muonSV:", total_pions_4muSV_2ph_2muSV)

# Debug category: events where we reconstruct two dimuonSVs but no 4muonSV
n_events_problematic = df.Filter("nPions_2ph_2muSV > 0 &&nPions_4muSV == 0").Count().GetValue()
print("Events with two matched dimuonSVs but NO 4muonSV:", n_events_problematic)
