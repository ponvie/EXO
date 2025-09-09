import numpy as np
import matplotlib.pyplot as plt

numHist, bin_edges = np.histogram(darkPhoton_pt_matched, bins=90)
denomHist, _ = np.histogram(darkPhoton_pt_gen, bins=90)

# Bin centers for plotting
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Efficiency
eff = []
for num, den in zip(numHist, denomHist):
    if den == 0:
        eff.append(0)
    else:
        val = num / den
        eff.append(val)

eff = np.array(eff)


plt.figure(figsize=(7,5))
plt.hist(darkPhoton_pt_gen, bins=90, histtype="step", linewidth=2, label="Generated dark photons")
plt.hist(darkPhoton_pt_matched, bins=90, histtype="step", linewidth=2, label="Matched dark photons")

plt.errorbar(bin_centers, eff, fmt='o', color="g", capsize=4, label="Efficiency")

plt.xlabel("Dark photon $p_T$ [GeV]")
plt.ylabel("Events / Efficiency")
plt.legend()
plt.savefig("Efficiency_darkPhoton_PT_RATIO.png")
plt.show()
