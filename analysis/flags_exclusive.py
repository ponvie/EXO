"""
if (event satisfies fourmuonSV selection):
     isFourMuonSV = true
else:
    if (event satisfies multimuonSV selection):
        isMultiMuonSV = true
    else:
        if (event satisfies singlemuonSV selection):
            isSingleMuonSV = true
"""


import ROOT
import os
import argparse

# --- Argomenti da condor (input e output) --------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input ROOT file")
parser.add_argument("--output_txt", required=True, help="Output TXT file")
parser.add_argument("--output_json", required=True, help="Output JSON file")
parser.add_argument("--output_root", required=True, help="Output ROOT file")
args = parser.parse_args()

# --- Carica file ROOT -----------------------------------------
df = ROOT.RDataFrame("Events", args.input)

# --- PDG IDs and parameters ------------------------------------------------
pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
deltaR_threshold = 0.05
muon_mass        = 0.105

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

// SV selection function
SVSelectionResult SelectSV(
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
    SVSelectionResult result{false, false, false, {}};

    // Step 1: Check four-muon SV with chi2 < 10
    for (size_t i = 0; i < fourmuonSV_chi2.size(); ++i) {
        if (fourmuonSV_chi2[i] < 10) {
            result.isFourMuonSV = true;
            result.muonIndices = {fourmuonSV_mu1index[i], fourmuonSV_mu2index[i],
                                  fourmuonSV_mu3index[i], fourmuonSV_mu4index[i]};
            return result;
        }
    }

    // Step 2: Filter two-muon SVs with chi2 < 10 and ΔR < 1.2
    struct SV {
        int idx;
        float chi2;
        float mass;
    };
    std::vector<SV> goodSVs;
    for (size_t i = 0; i < muonSV_chi2.size(); ++i) {
        if (muonSV_chi2[i] >= 10) continue;
        if (deltaR(muonSV_mu1eta[i], muonSV_mu1phi[i],
                   muonSV_mu2eta[i], muonSV_mu2phi[i]) >= 1.2) continue;

        goodSVs.push_back({(int)i, muonSV_chi2[i], muonSV_mass[i]});
    }

    // Step 2a: Try to pair SVs with mass difference < 3%
    int bestPair1 = -1, bestPair2 = -1;
    float bestChi2 = std::numeric_limits<float>::max();
    for (size_t i = 0; i < goodSVs.size(); ++i) {
        for (size_t j = i+1; j < goodSVs.size(); ++j) {
            float relDiff = std::fabs(goodSVs[i].mass - goodSVs[j].mass) / goodSVs[i].mass;
            if (relDiff < 0.03) {
                float minChi2 = std::min(goodSVs[i].chi2, goodSVs[j].chi2);
                if (minChi2 < bestChi2) {
                    bestChi2 = minChi2;
                    bestPair1 = goodSVs[i].idx;
                    bestPair2 = goodSVs[j].idx;
                }
            }
        }
    }

    if (bestPair1 != -1) {
        result.isMultiMuonSV = true;
        result.muonIndices = {
            muonSV_mu1index[bestPair1], muonSV_mu2index[bestPair1],
            muonSV_mu1index[bestPair2], muonSV_mu2index[bestPair2]};
        return result;
    }

    // Step 2b: If no pair, take the SV with smallest chi2
    if (!goodSVs.empty()) {
        int bestIdx = -1;
        float minChi2 = std::numeric_limits<float>::max();
        for (auto& sv : goodSVs) {
            if (sv.chi2 < minChi2) {
                minChi2 = sv.chi2;
                bestIdx = sv.idx;
            }
        }
        if (bestIdx != -1) {
            result.isSingleMuonSV = true;
            result.muonIndices = {muonSV_mu1index[bestIdx], muonSV_mu2index[bestIdx]};
            return result;
        }
    }

    return result; // nothing found
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

#nEvents = df.Count().GetValue()
#print("Total number of events:", nEvents)
#------------------------------------------------------------------------------------------------------#
#HLT general mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("MuonFiredHLT_mask",
               "MuonsFiredHLT(MuonBPark_fired_HLT_Mu10_Barrel_L1HP11_IP6)")

#------------------------------------------------------------------------------------------------------#
#Before loose mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("svResult", "SelectSV(fourmuonSV_chi2, fourmuonSV_mu1index, fourmuonSV_mu2index, fourmuonSV_mu3index, fourmuonSV_mu4index, muonSV_chi2, muonSV_mass, muonSV_mu1eta, muonSV_mu1phi, muonSV_mu2eta, muonSV_mu2phi, muonSV_mu1index, muonSV_mu2index)")

df = df.Define("isFourMuonSV", "svResult.isFourMuonSV")
df = df.Define("isMultiMuonSV", "svResult.isMultiMuonSV")
df = df.Define("isSingleMuonSV", "svResult.isSingleMuonSV")

total_fourmuonsSV     = df.Filter("isFourMuonSV").Count()
total_multimuonsSV    = df.Filter("isMultiMuonSV").Count()
total_singlemuonSV    = df.Filter("isSingleMuonSV").Count()
#------------------------------------------------------------------------------------------------------#
#After loose mask
#------------------------------------------------------------------------------------------------------#

df = df.Define("all_muons_are_loose",
               "SVAllLoose(MuonBPark_looseId, svResult.muonIndices)")

df = df.Define("isFourMuonSV_LOOSE", "(isFourMuonSV && all_muons_are_loose)")
df = df.Define("isMultiMuonSV_LOOSE", "(isMultiMuonSV && all_muons_are_loose)")
df = df.Define("isSingleMuonSV_LOOSE", "(isSingleMuonSV && all_muons_are_loose)")

loose_fourmuonsSV     = df.Filter("isFourMuonSV_LOOSE").Count()
loose_multimuonsSV    = df.Filter("isMultiMuonSV_LOOSE").Count()
loose_singlemuonSV    = df.Filter("isSingleMuonSV_LOOSE").Count()

#------------------------------------------------------------------------------------------------------#
#At least one muon fired HLT mask
#------------------------------------------------------------------------------------------------------#
df = df.Define("at_least_one_HLT",
               "AtLeastOneHLT(MuonFiredHLT_mask, svResult.muonIndices)")

df = df.Define("isFourMuonSV_LOOSE_HLT", "(isFourMuonSV && all_muons_are_loose && at_least_one_HLT)")
df = df.Define("isMultiMuonSV_LOOSE_HLT", "(isMultiMuonSV && all_muons_are_loose && at_least_one_HLT)")
df = df.Define("isSingleMuonSV_LOOSE_HLT", "(isSingleMuonSV && all_muons_are_loose && at_least_one_HLT)")

loose_HLT_fourmuonsSV   = df.Filter("isFourMuonSV_LOOSE_HLT").Count()
loose_HLT_multimuonsSV  = df.Filter("isMultiMuonSV_LOOSE_HLT").Count()
loose_HLT_singlemuonSV  = df.Filter("isSingleMuonSV_LOOSE_HLT").Count()

bkg_events = 99864988
sig_events = 1000000
# ------------------------------------------------------------------------------------------------------
# Tabella
# ------------------------------------------------------------------------------------------------------
def print_table_row(label, count):
    val = count.GetValue()
    return f"{label:<35} {val:>10}   {100.0*val/bkg_events:6.2f}%"

lines = []
lines.append("Selections summary (per-event counts):")
lines.append(f"Total Events: , {bkg_events} ")
#lines.append(f"FiredHLT Events: , {filtered_events.GetValue()} ")
lines.append("="*60)
lines.append(f"{'Selection':<35} {'Events':>10}   {'Fraction':>8}")
lines.append("="*60)
lines.append(print_table_row("isFourMuonSV", total_fourmuonsSV))
lines.append(print_table_row("isMultiMuonSV", total_multimuonsSV))
lines.append(print_table_row("isSingleMuonSV", total_singlemuonSV))
lines.append("-"*60)
lines.append(print_table_row("isFourMuonSV LOOSE", loose_fourmuonsSV))
lines.append(print_table_row("isMultiMuonSV LOOSE", loose_multimuonsSV))
lines.append(print_table_row("isSingleMuonSV LOOSE", loose_singlemuonSV))
lines.append("-"*60)
lines.append(print_table_row("isFourMuonSV LOOSE + ≥1 HLT", loose_HLT_fourmuonsSV))
lines.append(print_table_row("isMultiMuonSV LOOSE + ≥1 HLT", loose_HLT_multimuonsSV))
lines.append(print_table_row("isSingleMuonSV LOOSE + ≥1 HLT", loose_HLT_singlemuonSV))
lines.append("="*60)

# Scrivi su file txt
with open(args.output_txt, "w") as f:
    for l in lines:
        f.write(l + "\n")

# ------------------------------------------------------------------------------------------------------
# JSON - JavaScript Object Notation
# ------------------------------------------------------------------------------------------------------      
        
import json

results = {
    "TotalEvents": int(bkg_events),
    #"FiredHLTEvents": int(filtered_events.GetValue()),
    "Selections": {
        "isFourMuonSV": int(total_fourmuonsSV.GetValue()),
        "isMultiMuonSV": int(total_multimuonsSV.GetValue()),
        "isSingleMuonSV": int(total_singlemuonSV.GetValue()),
        "isFourMuonSV_LOOSE": int(loose_fourmuonsSV.GetValue()),
        "isMultiMuonSV_LOOSE": int(loose_multimuonsSV.GetValue()),
        "isSingleMuonSV_LOOSE": int(loose_singlemuonSV.GetValue()),
        "isFourMuonSV_LOOSE_HLT": int(loose_HLT_fourmuonsSV.GetValue()),
        "isMultiMuonSV_LOOSE_HLT": int(loose_HLT_multimuonsSV.GetValue()),
        "isSingleMuonSV_LOOSE_HLT": int(loose_HLT_singlemuonSV.GetValue())
    }
}

with open(args.output_json, "w") as f:
    json.dump(results, f, indent=2)
# ------------------------------------------------------------------------------------------------------
# Scrivi nuovo ROOT file con variabili
# ------------------------------------------------------------------------------------------------------
all_columns = df.GetColumnNames()

cols_to_save = [c for c in all_columns if c.startswith("muonSV_") or c.startswith("fourmuonSV_")]

# aggiungi anche i flag che avevi già
cols_to_save += [
    "isFourMuonSV", "isMultiMuonSV", "isSingleMuonSV",
    "isFourMuonSV_LOOSE", "isMultiMuonSV_LOOSE", "isSingleMuonSV_LOOSE",
    "isFourMuonSV_LOOSE_HLT", "isMultiMuonSV_LOOSE_HLT", "isSingleMuonSV_LOOSE_HLT",
    "nmuonSV", "nfourmuonSV"
]
df.Snapshot("Events", args.output_root, cols_to_save)

