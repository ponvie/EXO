#Efficiency

import numpy as np
import matplotlib.pyplot as plt
import ROOT, os
import math

# folder and file
dir = "/Users/veronicapilloni/Desktop/EXO/"
file = "nano.root"
events = ROOT.RDataFrame("Events", os.path.join(dir, file))

# PDG IDs
pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
pdgID_muon       = 13

ROOT.gInterpreter.Declare("""
    TLorentzVector CreateLorentzVector(double pt, double eta, double phi, double mass) {
        TLorentzVector v;
        v.SetPtEtaPhiM(pt, eta, phi, mass);
        return v;
    }
""")

# function to bring deltaPhi in [-pi, pi]
def Phi_mpi_pi(x):
    PI = 3.14159
    while x >= PI:
        x -= 2 * PI
    while x < -PI:
        x += 2 * PI
    return x

# function to compute deltaR using Phi_mpi_pi
def deltaR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = Phi_mpi_pi(phi1 - phi2)
    return math.sqrt(deta*deta + dphi*dphi)

def prepare_dataframe(df_raw):
    return (
        df_raw
        .Define("PDGid", "GenPart_pdgId")
        
        .Define("nDarkPions",    f"Sum(GenPart_pdgId == {pdgID_darkpion})")
        .Define("nDarkPhotons",  f"Sum(GenPart_pdgId == {pdgID_darkphoton})")

    )

# prepared dataframe
df = prepare_dataframe(events)

# extract columns
pdg      = df.AsNumpy(["PDGid"])["PDGid"]
nDarkPho = df.AsNumpy(["nDarkPhotons"])["nDarkPhotons"]
nGenPart = df.AsNumpy(["nGenPart"])["nGenPart"]
nmuonSV = df.AsNumpy(["nmuonSV"])["nmuonSV"]

GenPart_eta = df.AsNumpy(["GenPart_eta"])["GenPart_eta"]
GenPart_phi = df.AsNumpy(["GenPart_phi"])["GenPart_phi"]

muon1pt  = df.AsNumpy(["muonSV_mu1pt"])["muonSV_mu1pt"]
muon1eta = df.AsNumpy(["muonSV_mu1eta"])["muonSV_mu1eta"]
muon1phi = df.AsNumpy(["muonSV_mu1phi"])["muonSV_mu1phi"]

muon2pt  = df.AsNumpy(["muonSV_mu2pt"])["muonSV_mu2pt"]
muon2eta = df.AsNumpy(["muonSV_mu2eta"])["muonSV_mu2eta"]
muon2phi = df.AsNumpy(["muonSV_mu2phi"])["muonSV_mu2phi"]

print("dimension nGenPart:", len(nGenPart))
print("dimension nmuonSV:", len(nmuonSV))
print("\n")

deltaR_threshold = 0.5  
n_matched_photons = 0
for event in range(len(nGenPart)):
    
    mu1pt  = muon1pt[event]
    mu1eta = muon1eta[event]
    mu1phi = muon1phi[event]

    mu2pt  = muon2pt[event]
    mu2eta = muon2eta[event]
    mu2phi = muon2phi[event]

    dimuon_eta = []
    dimuon_phi = []
    for pair in range(len(mu1pt)):
        muon1 = ROOT.CreateLorentzVector(mu1pt[pair], mu1eta[pair], mu1phi[pair], 0.105)
        muon2 = ROOT.CreateLorentzVector(mu2pt[pair], mu2eta[pair], mu2phi[pair], 0.105)
        dimuon = muon1 + muon2
        dimuon_eta.append(dimuon.Eta())
        dimuon_phi.append(dimuon.Phi())

    #print(f"Event {event} dimuon_eta: {dimuon_eta}")
    #print(f"Event {event} dimuon_phi: {dimuon_phi}")
    
    for particle in range(nGenPart[event]):
        
        # skip particles different from dark photon
        if pdg[event][particle] != pdgID_darkphoton:
            continue
        
        darkphoton_eta = GenPart_eta[event][particle]
        darkphoton_phi = GenPart_phi[event][particle]
        
        # loop over reconstructed vertices and check deltaR with generated dark photons
        matched=False
        for vertex in range(len(dimuon_eta)):
            dR = deltaR(dimuon_eta[vertex], dimuon_phi[vertex], darkphoton_eta, darkphoton_phi)
            if dR < deltaR_threshold:
                matched=True
                break   
        if matched:
            n_matched_photons += 1
                           
print("Total generated dark photons:", nDarkPho.sum())
print("Number of matched dark photons:", n_matched_photons)

Efficiency = n_matched_photons / nDarkPho.sum() if nDarkPho.sum() > 0 else 0
Efficiency_percent = Efficiency * 100
print(f"Matching efficiency: {Efficiency:.4f}")
print(f"Matching efficiency (percent): {Efficiency_percent:.2f}%")
