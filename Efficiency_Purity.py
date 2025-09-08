import numpy as np
import matplotlib.pyplot as plt
import ROOT, os

dir = "/Users/veronicapilloni/Desktop/EXO/"
file ="nano.root"
events = ROOT.RDataFrame("Events", os.path.join(dir, file))


pdgID_darkphoton = 999999
pdgID_darkpion   = 4900111
pdgID_muon       = 13


def prepare_dataframe(df_raw):
    return (
        df_raw
        .Define("gen_idx", "nGenPart")  #Index into genParticle list for MC matching to status==1 muons
        .Define("vertex_idx", "Muon_genPartIdx") #Index into genParticle list for MC matching to status==1 muons
        .Define("PDGid", "GenPart_pdgId") #PDG ID of the particle
        .Define("MotherIdx", "GenPart_genPartIdxMother") #PDG ID of the mother particle
        
        .Define("nDarkPions", f"Sum(GenPart_pdgId == {pdgID_darkpion})") #tot gen dark pions
        .Define("nDarkPhotons", f"Sum(GenPart_pdgId == {pdgID_darkphoton})") #tot gen dark photons
        
    )

df = prepare_dataframe(events)

vertex = df.AsNumpy(["vertex_idx"])["vertex_idx"]
gen  = df.AsNumpy(["gen_idx"])["gen_idx"]
pdg  = df.AsNumpy(["PDGid"])["PDGid"]
mother = df.AsNumpy(["MotherIdx"])["MotherIdx"]

print("dimensione gen_idx:", len(gen))
print("dimensione vertex_idx:", len(vertex))
print("dimensione PDGid:", len(pdg))
print("dimensione MotherIdx:", len(mother))
print("\n")

print("pdgID[0]:", pdg[0][125])
print("motherID[0]:", mother[0])
for i in range(5):
    print(f"gen_idx[{i}]:", gen[i])
    print(f"vertex_idx[{i}]:", vertex[i])
    print("\n")

two_muon_vertices = 0
four_muon_vertices = 0


# iterate over reconstructed vertices
for v in range(len(vertex)):
    n_eff_muons = 0
    print(f"vertex_idx[{v}]:", vertex[v])
    # keep only those with at least two "muons"
    if len(vertex[v]) < 2:
        print("less than 2 muons")
    
        continue
    # iterate over the particles associated with the vertex
    for v_muon in range(len(vertex[v])):
        if abs(pdg[v][v_muon]) != pdgID_muon:
            print(len(vertex[v]))
            print("not a muon")
            continue

        mother_idx = mother[v][v_muon]
        if mother_idx < 0:
            continue
        
        #only consider muons coming from dark photons
        if pdg[v][mother_idx] == pdgID_darkphoton:
            n_eff_muons += 1
    
    if n_eff_muons == 2:
        two_muon_vertices += 1
        
    elif n_eff_muons == 4:
        four_muon_vertices += 1


print("Number of 2-muon vertices:", two_muon_vertices)
print("Number of 4-muon vertices:", four_muon_vertices)
