"""
if (event satisfies fourmuonSV selection):
     isFourMuonSV = true
if (event satisfies  multimuonSV selection):
    isMultiMuonSV = true
if (event satisfies  singlemuonSV selection):
    isSingleMuonSV = true
"""


import ROOT, os

# --- Input file ------------------------------------------------------------
dir  = "/Users/veronicapilloni/Desktop/EXO/newfile/"
file = "data_0.root"
df = ROOT.RDataFrame("Events", os.path.join(dir, file))

# --- PDG IDs and parameters ------------------------------------------------
pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
deltaR_threshold = 0.5
muon_mass        = 0.105

nEvents = df.Count().GetValue()
print("Total number of events:", nEvents)


# Step 1: Select events where the HLT path fired
# Mu10   → requires a muon with pT > 10 GeV
# Barrel → restricts the muon to the barrel region of the detector, i.e. |η| < ~1.0–1.2
# L1HP11 → the L1 seed used: Level-1 High-Purity single-muon with pT threshold ~11 GeV
# IP6    → means the muon must have impact parameter > 0.06 cm (i.e. 600 µm)

#df_HLT = df.Filter("HLT_Mu10_Barrel_L1HP11_IP6")
#nEvents_fired = df_HLT.Count().GetValue()
#print("Total number of events after HLT:", nEvents_fired)

#OR

ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
#include <cmath>
#include <limits>
using ROOT::VecOps::RVec;

// Function that returns a mask (RVec<bool>) telling which muons fired the HLT
RVec<bool> MuonsFiredHLT(const RVec<int>& firedHLT) {
    RVec<bool> mask(firedHLT.size(), false);
    for (size_t i = 0; i < firedHLT.size(); ++i) {
        if (firedHLT[i] == 1) {
            mask[i] = true;
        }
    }
    return mask;
    }

struct SVSelectionResult {
    bool isFourMuonSV;
    bool isMultiMuonSV;
    bool isSingleMuonSV;

    int selectedFourMuonSV;  // index of chosen four-muon SV (-1 if none)
    std::pair<int,int> selectedTwoMuonSVs; // indices of the chosen pair (-1,-1 if none)
    int selectedSingleMuonSV; // index of chosen single SV (-1 if none)

    RVec<int> muonIndices; // optional: keep all muon indices used
};

// Check ΔR
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2*M_PI - dphi;
    return std::sqrt(deta*deta + dphi*dphi);
}

SVSelectionResult SelectSV_independent(
    const RVec<float>& fourmuonSV_chi2,
    const RVec<int>& fourmuonSV_mu1index,
    const RVec<int>& fourmuonSV_mu2index,
    const RVec<int>& fourmuonSV_mu3index,
    const RVec<int>& fourmuonSV_mu4index,
    const RVec<float>& muonSV_chi2,
    const RVec<float>& muonSV_mass,
    const RVec<float>& muonSV_mu1eta,
    const RVec<float>& muonSV_mu1phi,
    const RVec<float>& muonSV_mu2eta,
    const RVec<float>& muonSV_mu2phi,
    const RVec<int>& muonSV_mu1index,
    const RVec<int>& muonSV_mu2index
) {
    SVSelectionResult result{false, false, false, -1, {-1,-1}, -1, {}};

    // --- Four-muon SV check ---
    for (size_t i = 0; i < fourmuonSV_chi2.size(); ++i) {
        if (fourmuonSV_chi2[i] < 10) {
            result.isFourMuonSV = true;
            result.selectedFourMuonSV = i;
            result.muonIndices.insert(result.muonIndices.end(),
                {fourmuonSV_mu1index[i], fourmuonSV_mu2index[i],
                 fourmuonSV_mu3index[i], fourmuonSV_mu4index[i]});
        }
    }

    // --- Two-muon SV selection ---
    struct SV { int idx; float chi2; float mass; };
    std::vector<SV> goodSVs;
    for (size_t i = 0; i < muonSV_chi2.size(); ++i) {
        if (muonSV_chi2[i] >= 10) continue;
        if (deltaR(muonSV_mu1eta[i], muonSV_mu1phi[i],
                   muonSV_mu2eta[i], muonSV_mu2phi[i]) >= 1.2) continue;
        goodSVs.push_back({(int)i, muonSV_chi2[i], muonSV_mass[i]});
    }

    // Try to pair SVs with mass difference < 3%
    for (size_t i = 0; i < goodSVs.size(); ++i) {
        for (size_t j = i+1; j < goodSVs.size(); ++j) {
            float relDiff = std::fabs(goodSVs[i].mass - goodSVs[j].mass) / goodSVs[i].mass;
            if (relDiff < 0.03) {
                result.isMultiMuonSV = true;
                result.selectedTwoMuonSVs = {goodSVs[i].idx, goodSVs[j].idx};
                result.muonIndices.insert(result.muonIndices.end(),
                    {muonSV_mu1index[goodSVs[i].idx], muonSV_mu2index[goodSVs[i].idx],
                     muonSV_mu1index[goodSVs[j].idx], muonSV_mu2index[goodSVs[j].idx]});
            }
        }
    }

    // Single SV: best chi2
    if (!goodSVs.empty()) {
        result.isSingleMuonSV = true;
        auto bestIt = std::min_element(
            goodSVs.begin(), goodSVs.end(),
            [](const SV& a, const SV& b){ return a.chi2 < b.chi2; });
        result.selectedSingleMuonSV = bestIt->idx;
        result.muonIndices.insert(result.muonIndices.end(),
            {muonSV_mu1index[bestIt->idx], muonSV_mu2index[bestIt->idx]});
    }

    return result;
}

// Generic function: all muons must be Loose and fire HLT
bool SVAllLoose(
    const RVec<int>& muonLooseId,
    const RVec<int>& muonIdxs
) {
    for (int idx : muonIdxs) {
        if (idx < 0 || idx >= (int)muonLooseId.size()) return false;
        if (muonLooseId[idx] != 1) return false; // must be Loose
    }
    return true;
}

// At least one muon from the selected SV must fire HLT
bool AtLeastOneHLT(
    const RVec<bool>& muonHLTmask,
    const RVec<int>& muonIdxs
) {
    for (int idx : muonIdxs) {
        if (idx < 0 || idx >= (int)muonHLTmask.size()) return false;
        if (muonHLTmask[idx]) return true;  // found at least one
    }
    return false;
}


""")

#------------------------------------------------------------------------------------------------------#
#HLT general mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("MuonFiredHLT_mask",
               "MuonsFiredHLT(MuonBPark_fired_HLT_Mu10_Barrel_L1HP11_IP6)")

df = df.Define("nMuons", "Muon_pt.size()")  
total_muons = df.Sum("nMuons")

df = df.Define("nSelectedMuons", "Muon_pt[MuonFiredHLT_mask].size()")
selected_muons = df.Sum("nSelectedMuons")



#------------------------------------------------------------------------------------------------------#
#Before loose mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("svResult", "SelectSV_independent(fourmuonSV_chi2, fourmuonSV_mu1index, fourmuonSV_mu2index, fourmuonSV_mu3index, fourmuonSV_mu4index, muonSV_chi2, muonSV_mass, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2eta, muonSV_mu2phi, muonSV_mu1index, muonSV_mu2index)")

df = df.Define("isFourMuonSV", "svResult.isFourMuonSV")
total_fourmuonsSV = df.Sum("isFourMuonSV")

df = df.Define("isMultiMuonSV", "svResult.isMultiMuonSV")
total_multimuonsSV = df.Sum("isMultiMuonSV")

df = df.Define("isSingleMuonSV", "svResult.isSingleMuonSV")
total_singlemuonSV = df.Sum("isSingleMuonSV")



#------------------------------------------------------------------------------------------------------#
#After loose mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("all_muons_are_loose",
               "SVAllLoose(MuonBPark_looseId, svResult.muonIndices)")

df = df.Define("isFourMuonSV_LOOSE", "(isFourMuonSV && all_muons_are_loose)")
loose_fourmuonsSV = df.Sum("isFourMuonSV_LOOSE")

df = df.Define("isMultiMuonSV_LOOSE", "(isMultiMuonSV && all_muons_are_loose)")
loose_multimuonsSV = df.Sum("isMultiMuonSV_LOOSE")

df = df.Define("isSingleMuonSV_LOOSE", "(isSingleMuonSV && all_muons_are_loose)")
loose_singlemuonSV = df.Sum("isSingleMuonSV_LOOSE")


#------------------------------------------------------------------------------------------------------#
#At least one muon fired HLT mask
#------------------------------------------------------------------------------------------------------#
df = df.Define("at_least_one_HLT",
               "AtLeastOneHLT(MuonFiredHLT_mask, svResult.muonIndices)")

df = df.Define("isFourMuonSV_LOOSE_HLT", "(isFourMuonSV && all_muons_are_loose && at_least_one_HLT)")
loose_HLT_fourmuonsSV = df.Sum("isFourMuonSV_LOOSE")


df = df.Define("isMultiMuonSV_LOOSE_HLT", "(isMultiMuonSV && all_muons_are_loose && at_least_one_HLT)")
loose_HLT_multimuonsSV = df.Sum("isMultiMuonSV_LOOSE_HLT")


df = df.Define("isSingleMuonSV_LOOSE_HLT", "(isSingleMuonSV && all_muons_are_loose && at_least_one_HLT)")
loose_HLT_singlemuonSV = df.Sum("isSingleMuonSV_LOOSE_HLT")



print("Total muons before mask:", total_muons.GetValue())
print("Total muons after mask:", selected_muons.GetValue())
print("\n")
print("isFourMuonSV:", total_fourmuonsSV.GetValue())
print("isMultiMuonSV:", total_multimuonsSV.GetValue())
print("isSingleMuonSV:", total_singlemuonSV.GetValue())
print("\n")
print("isFourMuonSV LOOSE:", loose_fourmuonsSV.GetValue())
print("isMultiMuonSV LOOSE:", loose_multimuonsSV.GetValue())
print("isSingleMuonSV LOOSE:", loose_singlemuonSV.GetValue())
print("\n")
print("isFourMuonSV LOOSE + 1 muon fired HLT:", loose_HLT_fourmuonsSV.GetValue())
print("isMultiMuonSV LOOSE + 1 muon fired HLT:", loose_HLT_multimuonsSV.GetValue())
print("isSingleMuonSV LOOSE + 1 muon fired HLT:", loose_HLT_singlemuonSV.GetValue())
