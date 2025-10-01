# Run after ROOT_MuonSVs.py
# --- Declare C++ helper functions ---
ROOT.gInterpreter.Declare("""

#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <cmath>
using namespace ROOT::VecOps;

// DeltaR helper
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    float deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// PDG ID → carica
int pdgToCharge(int pdgId) {
    if (pdgId == 13)  return -1;  // mu-
    if (pdgId == -13) return +1;  // mu+
    return 0;
}

// --- Funzione 1: flag = 1 se muone viene da un dark photon ---
RVec<int> PhotonMuonMask(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    int pdgID_darkphoton
) {
    RVec<int> mask(pdgId.size(), 0);

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) != 13) continue; // solo muoni

        int mom = motherIdx[i];
        if (mom < 0 || mom >= (int)pdgId.size()) continue;

        if (pdgId[mom] == pdgID_darkphoton) {
            mask[i] = 1;
        }
    }
    return mask;
}

// --- Funzione 2: flag = 1 se carica del reco muone matcha quella del gen muone ---
RVec<int> ChargeMatchMask(
    const RVec<int>& pdgId,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    const RVec<float>& reco_eta,
    const RVec<float>& reco_phi,
    const RVec<int>& reco_charge,
    float dr_threshold
) {
    RVec<int> mask(pdgId.size(), 0);

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) != 13) continue;

        int true_charge = pdgToCharge(pdgId[i]);
        float best_dr = 999.;
        int best_idx = -1;

        for (size_t j = 0; j < reco_eta.size(); ++j) {
            float dr = deltaR(gp_eta[i], gp_phi[i], reco_eta[j], reco_phi[j]);
            if (dr < best_dr && dr < dr_threshold) {
                best_dr = dr;
                best_idx = j;
            }
        }

        if (best_idx >= 0 && reco_charge[best_idx] == true_charge) {
            mask[i] = 1;
        }
    }
    return mask;
}
""")

df = (
    df
    # Muoni provenienti da dark photon
    .Define("mask_matchedPhotonMuon",
            f"PhotonMuonMask(GenPart_pdgId, GenPart_genPartIdxMother, {pdgID_darkphoton})")

    # Muoni con carica corretta rispetto al gen
    .Define("mask_chargeMatched",
            f"ChargeMatchMask(GenPart_pdgId, GenPart_eta, GenPart_phi, Muon_eta, Muon_phi, Muon_charge, {deltaR_threshold})")

    # Muoni che soddisfano entrambe le condizioni (esplicito ==1)
    .Define("mask_goodMuons", "(mask_matchedPhotonMuon == 1) && (mask_chargeMatched == 1)")

    # conteggi
    .Define("n_muons_fromPhoton", "Sum(mask_matchedPhotonMuon==1)")
    .Define("n_goodMuons", "Sum(mask_goodMuons)")
)

print("Total muons from photons:", df.Sum("n_muons_fromPhoton").GetValue())
print("Good muons (from photon & correct charge):", df.Sum("n_goodMuons").GetValue())
