numHist = np.histogram(darkPhoton_pt_matched, bins=90)[0]
denomHist  = np.histogram(darkPhoton_pt_gen, bins=90)[0]

print(numHist)
print(denomHist)

eff=[]
for num, den in zip(numHist, denomHist):
    if den ==0:
        eff.append(0)
    else:
        eff.append(num/den)
        
plt.figure(figsize=(7,5))
plt.hist(darkPhoton_pt_gen, bins=80, histtype="step", linewidth=2, label="Generated dark photons")
plt.hist(darkPhoton_pt_matched, bins=80, histtype="step", linewidth=2, label="Matched dark photons")
plt.hist(eff, bins=80, histtype="step", linewidth=2, label="Ratio", color="g")
plt.xlabel("Dark photon $p_T$ [GeV]")
plt.ylabel("Events")
plt.xlabel("Dark photon $p_T$ [GeV]")
plt.ylabel("Events")
plt.legend()
plt.savefig("Efficiency_darkPhoton_PT_RATIO.png")
plt.show()
