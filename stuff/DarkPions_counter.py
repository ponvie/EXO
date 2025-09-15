import ROOT, os
import matplotlib.pyplot as plt

dir = "/Users/veronicapilloni/Desktop/EXO/"
file = "data_1.root"

f = ROOT.TFile(os.path.join(dir, file))
tree = f.Get("Events")

pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
pdgID_muon       = 13

nPions_2mu = []
nPions_4mu = []

for event in tree:
    nPart = event.nGenPart
    
    
    pion_muons = {}

    for i in range(nPart):
        if abs(event.GenPart_pdgId[i]) != pdgID_muon:
            continue

        
        mother_idx = event.GenPart_genPartIdxMother[i]
        if mother_idx < 0: 
            continue
        if event.GenPart_pdgId[mother_idx] != pdgID_darkphoton:
            continue

        pion_idx = event.GenPart_genPartIdxMother[mother_idx]
        if pion_idx < 0: 
            continue
        if event.GenPart_pdgId[pion_idx] != pdgID_darkpion:
            continue

        if pion_idx not in pion_muons:
            pion_muons[pion_idx] = 0
        pion_muons[pion_idx] += 1

    pions_2mu = sum(1 for nmu in pion_muons.values() if nmu >= 2)
    pions_4mu = sum(1 for nmu in pion_muons.values() if nmu == 4)

    nPions_2mu.append(pions_2mu)
    nPions_4mu.append(pions_4mu)

plt.hist(nPions_2mu, bins=20, histtype="step", linewidth=2, label="Dark pions → 2μ")
plt.hist(nPions_4mu, bins=20, histtype="step", linewidth=2, label="Dark pions → 4μ")
plt.xlabel("Number of dark pions")
plt.ylabel("Events")
plt.legend()
plt.show()
