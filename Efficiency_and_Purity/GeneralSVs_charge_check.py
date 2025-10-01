# Run after ROOT_muonSVs.py
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <cmath>
using namespace ROOT::VecOps;

// DeltaR helper
float deltaR(float eta1, float phi1, float eta2, float phi2){
    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    float deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// PDG → charge
int pdgToCharge(int pdgId){
    if(pdgId==13) return -1;
    if(pdgId==-13) return +1;
    return 0;
}


ROOT::RVec<int> MatchRecoMuons(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& motherIdx,
    const ROOT::RVec<float>& gp_eta,
    const ROOT::RVec<float>& gp_phi,
    const ROOT::RVec<float>& reco_eta,
    const ROOT::RVec<float>& reco_phi,
    int targetMotherPDG,
    bool requireGrandmother=false,
    int grandmotherPDG=0,
    float dr_threshold=0.03
){
    ROOT::RVec<int> mask(reco_eta.size(),0);

    for(size_t i=0;i<pdgId.size();++i){
        if(abs(pdgId[i])!=13) continue; // solo muoni

        int mom = motherIdx[i];
        if(mom<0 || mom>=(int)pdgId.size()) continue;
        if(pdgId[mom]!=targetMotherPDG) continue; // madre = dark photon

        // controllo sul nonno se richiesto
        if(requireGrandmother){
            int grandma = motherIdx[mom]; // madre del mother
            if(grandma<0 || grandma>=(int)pdgId.size()) continue;
            if(pdgId[grandma]!=grandmotherPDG) continue;
        }

        // match al reco
        float best_dr = 999.;
        int best_idx = -1;
        for(size_t j=0;j<reco_eta.size();++j){
            float dr = deltaR(gp_eta[i],gp_phi[i],reco_eta[j],reco_phi[j]);
            if(dr<best_dr && dr<dr_threshold){
                best_dr = dr;
                best_idx = j;
            }
        }
        if(best_idx>=0) mask[best_idx]=1;
    }
    return mask;
}

// Charge matching condizionata sui reco già matched
ROOT::RVec<int> ChargeMatchedReco(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& motherIdx,
    const ROOT::RVec<float>& gp_eta,
    const ROOT::RVec<float>& gp_phi,
    const ROOT::RVec<float>& reco_eta,
    const ROOT::RVec<float>& reco_phi,
    const ROOT::RVec<int>& reco_charge,
    const ROOT::RVec<int>& matched_mask,
    int targetMotherPDG,
    bool requireGrandmother=false,
    int grandmotherPDG=0,
    float dr_threshold=0.03
){
    ROOT::RVec<int> mask(reco_eta.size(),0);

    for(size_t i=0;i<pdgId.size();++i){
        if(abs(pdgId[i])!=13) continue;

        int mom = motherIdx[i];
        if(mom<0 || mom>=(int)pdgId.size()) continue;
        if(pdgId[mom]!=targetMotherPDG) continue;

        if(requireGrandmother){
            int grandma = motherIdx[mom];
            if(grandma<0 || grandma>=(int)pdgId.size()) continue;
            if(pdgId[grandma]!=grandmotherPDG) continue;
        }

        int true_charge = pdgToCharge(pdgId[i]);

        // match condizionato sui già matched
        float best_dr=999.;
        int best_idx=-1;
        for(size_t j=0;j<reco_eta.size();++j){
            if(matched_mask[j]==0) continue;
            float dr = deltaR(gp_eta[i],gp_phi[i],reco_eta[j],reco_phi[j]);
            if(dr<best_dr && dr<dr_threshold){
                best_dr=dr;
                best_idx=j;
            }
        }
        if(best_idx>=0 && reco_charge[best_idx]==true_charge)
            mask[best_idx]=1;
    }

    return mask;
}
""")

df = (
    df
    # Muoni da dark photon
    .Define("mask_reco_fromDarkPhoton",
            f"MatchRecoMuons(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, Muon_eta, Muon_phi, {pdgID_darkphoton}, false, 0, {deltaR_threshold})")
    .Define("mask_chargeMatched_darkphoton",
            f"ChargeMatchedReco(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, Muon_eta, Muon_phi, Muon_charge, mask_reco_fromDarkPhoton, {pdgID_darkphoton}, false, 0, {deltaR_threshold})")
    .Define("n_reco_fromDarkPhoton", "Sum(mask_reco_fromDarkPhoton)")
    .Define("n_reco_good_darkphoton", "Sum(mask_chargeMatched_darkphoton)")

    # Muoni da dark pion -> dark photon
    .Define("mask_reco_fromPionPhoton",
            f"MatchRecoMuons(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, Muon_eta, Muon_phi, {pdgID_darkphoton}, true, {pdgID_darkpion}, {deltaR_threshold})")
    .Define("mask_chargeMatched_PionPhoton",
            f"ChargeMatchedReco(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, Muon_eta, Muon_phi, Muon_charge, mask_reco_fromPionPhoton, {pdgID_darkphoton}, true, {pdgID_darkpion}, {deltaR_threshold})")
    .Define("n_reco_fromPionPhoton", "Sum(mask_reco_fromPionPhoton==1)")
    .Define("n_reco_good_PionPhoton", "Sum(mask_chargeMatched_PionPhoton==1)")
)


n_reco_fromDarkPhoton     = df.Sum("n_reco_fromDarkPhoton").GetValue()
n_reco_good_darkphoton    = df.Sum("n_reco_good_darkphoton").GetValue()
n_reco_fromPionPhoton     = df.Sum("n_reco_fromPionPhoton").GetValue()
n_reco_good_PionPhoton    = df.Sum("n_reco_good_PionPhoton").GetValue()

print(f"Reco muons from dark photon: {n_reco_fromDarkPhoton}")
print(f"Reco muons with correct charge (dark photon): {n_reco_good_darkphoton}")
print(f"Reco muons from pion->photon: {n_reco_fromPionPhoton}")
print(f"Reco muons with correct charge (pion->photon): {n_reco_good_PionPhoton}")
