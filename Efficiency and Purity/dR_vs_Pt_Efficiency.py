# add to ROOT_muonSV.py before running 

import ROOT
from plotting_tools.root import get_labels, Canvas

ROOT.gStyle.SetOptStat(0)

# --- C++: Δphi e funzione che ritorna un vettore allineato a GenPart ---
ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
#include "ROOT/RVec.hxx"
using ROOT::VecOps::RVec;

static inline double deltaPhi_safe(double a, double b) {
    double d = a - b;
    while (d > M_PI)  d -= 2*M_PI;
    while (d <= -M_PI) d += 2*M_PI;
    return d;
}

// Returns a vector with size == pdgId.size(), drs[ip] = ΔR between the two muons
// daughters of gen particle ip if ip is a dark photon that decays to 2 muons,
// otherwise drs[ip] = -1
RVec<float> ComputeDeltaR_DarkPhotonMuons_perGenPart(
    const RVec<int>& pdgId,
    const RVec<int>& motherIdx,
    const RVec<float>& gp_eta,
    const RVec<float>& gp_phi,
    int pdgID_darkphoton,
    int pdgID_muon
) {
    size_t N = pdgId.size();
    RVec<float> drs(N, -1.0f);

    for (size_t ip = 0; ip < N; ++ip) {
        if (pdgId[ip] != pdgID_darkphoton) continue;

        int first = -1, second = -1;
        // find two muon daughters with mother == ip
        for (size_t j = 0; j < N; ++j) {
            if (std::abs(pdgId[j]) == pdgID_muon && motherIdx[j] == (int)ip) {
                if (first == -1) first = j;
                else { second = j; break; }
            }
        }
        if (first != -1 && second != -1) {
            double deta = double(gp_eta[first]) - double(gp_eta[second]);
            double dphi = deltaPhi_safe(double(gp_phi[first]), double(gp_phi[second]));
            drs[ip] = float(std::sqrt(deta*deta + dphi*dphi));
        }
    }

    return drs;
}
""")

# ---------------------- Python / RDF ----------------------
# assume df è il tuo RDataFrame già preparato e che esistano:
# - isDarkPhoton_into2mu (mask su GenPart con i dark photon che hanno 2 muoni)
# - darkPhoton_flags (mask su GenPart con i dark photon matched)
# - GenPart_* arrays (pdgId, genPartIdxMother, eta, phi)

pdgID_darkphoton = 999999
pdgID_muon = 13

# Define ΔR aligned to GenPart
df = df.Define(
    "DeltaR_perGenPart",
    f"ComputeDeltaR_DarkPhotonMuons_perGenPart(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, {pdgID_darkphoton}, {pdgID_muon})"
)

# Select only dark photons that actually have 2 muons (mask isDarkPhoton_into2mu must exist)
df = df.Define("DeltaR_darkPhoton_all", "DeltaR_perGenPart[isDarkPhoton_into2mu==1]")

# Select only those dark photons that are also matched (darkPhoton_flags has GenPart length)
df = df.Define("DeltaR_darkPhoton_matched", "DeltaR_perGenPart[(isDarkPhoton_into2mu==1) && (darkPhoton_flags==1)]")

# --- Istogrammi ---
h_all = df.Histo1D(("h_all", ";#Delta R(#mu,#mu) from dark photon;Events", 90, 0, 1), "DeltaR_darkPhoton_all").Clone()
h_matched = df.Histo1D(("h_matched", ";#Delta R(#mu,#mu) from dark photon;Events", 90, 0, 1), "DeltaR_darkPhoton_matched").Clone()

h_all.SetLineColor(ROOT.kBlue)
h_all.SetLineWidth(2)
h_matched.SetLineColor(ROOT.kRed)
h_matched.SetLineWidth(2)

# --- Plot CMS style ---
c = Canvas("c_dR", 1200, 700)
h_all.Draw("HIST")
h_matched.Draw("HIST SAME")

leg = ROOT.TLegend(0.60, 0.68, 0.90, 0.88)
leg.SetBorderSize(0)
leg.SetTextSize(0.024)
leg.AddEntry(h_all, "All dark photons", "l")
leg.AddEntry(h_matched, "Matched dark photons", "l")
leg.Draw()

# Labels CMS
draw_labels = get_labels(
    upper_left="Private Work",
    upper_right="2024 Simulation (13.6 TeV)",
    upper_left_offset=0.05
)
for label in draw_labels:
    label.Draw("same")

c.Update()
c.SaveAs("/Users/veronicapilloni/Desktop/EXO/data/darkPhoton_muonpair_deltaR.png")

print("Entries all (from hist):", h_all.Integral())
print("Entries matched (from hist):", h_matched.Integral())

ROOT.gStyle.SetOptStat(0)  # niente box statistico
# Istogrammi già creati: h_all, h_matched
eff = ROOT.TEfficiency(h_matched, h_all)

# Canvas
c2 = Canvas()

# Disegna efficienza
eff.SetLineColor(ROOT.kRed)
eff.Draw("AP")  # Axis + Points
eff.SetTitle(";#Delta R;Efficiency")
# Labels CMS
draw_labels = get_labels(
    upper_left="Private Work",
    upper_right="2024 Simulation (13.6 TeV)",
    upper_left_offset=0.05
)
for label in draw_labels:
    label.Draw("same")

# Aggiorna e salva
c2.Update()
c2.SaveAs("/Users/veronicapilloni/Desktop/EXO/data/EFF_DarkPhotons_dR_Gen_vs_Matched.png")


ROOT.gStyle.SetOptStat(0)  # disabilita box statistico

# Istogrammi 2D: All vs Matched

h2_all = df.Histo2D(
    ("h2_all", ";p_{T} [GeV];d_{xy} [cm]", 61, 0, 60, 50, 0, 1),
    "darkPhoton_pt", "DeltaR_darkPhoton_all"
).Clone()

h2_matched = df.Histo2D(
    ("h2_matched", ";p_{T} [GeV];d_{xy} [cm]", 61, 0, 60, 50, 0, 1),
    "darkPhoton_pt_matched", "DeltaR_darkPhoton_matched"
).Clone()


# Efficienza 2D
eff2D = ROOT.TEfficiency(h2_matched, h2_all)

# Canvas CMS style
c2 = Canvas()
c2.canvas.SetLeftMargin(0.16)
c2.canvas.SetRightMargin(0.20)
c2.canvas.SetBottomMargin(0.14)
eff2D.SetTitle(";p_{T} [GeV];#Delta R (#mu#mu);Efficiency")
eff2D.Draw("COLZ")  

# CMS labels
draw_labels = get_labels(
    upper_left="Private Work",
    upper_right="2024 Simulation (13.6 TeV)",
    upper_left_offset=0.05
)
for label in draw_labels:
    label.Draw("same")

c2.Update()
c2.SaveAs("/Users/veronicapilloni/Desktop/EXO/data/EFF_DarkPhotons_pt_vs_dR.png")

