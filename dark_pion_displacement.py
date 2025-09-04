import ROOT
import matplotlib.pyplot as plt
import os
import pandas as pd
from math import nan, isnan
import glob
from numpy import trapz
import numpy as np
import mplhep

dir = "/Users/veronicapilloni/Desktop/EXO/"
file ="data_1.root"

events = ROOT.RDataFrame("Events", os.path.join(dir, file))

pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
pdgID_muon       = 13

df = df.Filter(f"Sum(abs(GenPart_pdgId) == {pdgID_muon}) == 4")  # solo eventi con 4 muoni

#[particella, particella madre]
ID_arrays = df.AsNumpy(["GenPart_pdgId","GenPart_genPartIdxMother"])

#numero di muoni provenienti da dark photon per evento
nDPdecayToMu = []
for pdgIds, mothers in zip(ID_arrays["GenPart_pdgId"], ID_arrays["GenPart_genPartIdxMother"]):
    count = 0
    for i, pid in enumerate(pdgIds):
        if abs(pid) == 13:  # muone
            mom_idx = mothers[i]
            if mom_idx >= 0 and pdgIds[mom_idx] == pdgID_darkphoton: # proveniente da dark photon
                count += 1
    nDPdecayToMu.append(count)

# istogramma matplotlib
import matplotlib.pyplot as plt
plt.hist(nDPdecayToMu, bins=range(max(nDPdecayToMu)+2), histtype="step")
plt.xlabel("Number of muons from dark pions")
plt.ylabel("Events")
plt.show()
