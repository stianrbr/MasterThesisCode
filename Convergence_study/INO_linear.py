from Read_WADAM.misc import WADAM_res

import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


mainfolder = "S:\\Master\\Convergence study - body mesh\\"

vessel = "INO"

if vessel == "INO":
    title = "Peripheral tower design"
elif vessel == "Volturn":
    title = "Central tower design"
else:
    raise Exception("Check vessel")

mesh_0_8 = mainfolder+ vessel+"\\Mesh_0_8m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_0_8 = WADAM_res(filename=mesh_0_8)

mesh_1_2 = mainfolder+ vessel + "\\Mesh_1_2m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_2 = WADAM_res(filename=mesh_1_2)

mesh_1_8 = mainfolder+ vessel + "\\Mesh_1_8m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_8 = WADAM_res(filename=mesh_1_8)

mesh_2_7 = mainfolder+ vessel + "\\Mesh_2_7m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_2_7 = WADAM_res(filename=mesh_2_7)

results = [res_0_8, res_1_2, res_1_8, res_2_7]
labels = ["0.8 m", "1.2 m", "1.8 m", "2.7 m"]

try:
    os.mkdir("Comparison/INO")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/INO\\Added_mass")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/INO\\Damping")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/INO\\Motion_RAO")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/INO\\Load_RAO")
except FileExistsError:
    pass



for i in range(6):
    for j in range(6):

        added_mass, period = zip(*[res.get_Total_added_mass(i,j) for res in results])
        fig = plt.figure()
        ax = plt.axes()
        for a,p ,l in zip(added_mass, period, labels): ax.plot(p, a, label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Added mass - (i,j)=({i+1},{j+1})")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Added_mass/a_{i+1}{j+1}.png")
        plt.close(fig)



        damping, period = zip(*[res.get_Total_damping(i,j) for res in results])

        fig = plt.figure()
        ax = plt.axes()
        for b, p, l in zip(damping, period, labels): ax.plot(p, b, label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Damping - (i,j)=({i + 1},{j + 1})")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Damping/b_{i + 1}{j + 1}.png")
        plt.close(fig)

headings = results[0].Environment_data.wave_headings

for j, h in enumerate(headings):
    for i in range(6):
        amp, phase, period = zip(*[res.get_Motion_RAO(i,j) for res in results])

        fig = plt.figure()
        ax = plt.axes()
        for a, per, l in zip(amp, period, labels): ax.plot(per, a, label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Motion_RAO/x_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        for a, per, l in zip(phase, period, labels): ax.plot(per, a, label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Motion_RAO/x_ph_{i + 1}_head_{h}.png")
        plt.close(fig)

        amp, phase, period = zip(*[res.get_Load_RAO(i,j) for res in results])


        fig = plt.figure()
        ax = plt.axes()
        for a, per, l in zip(amp, period, labels): ax.plot(per, a, label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Load_RAO/f_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        for a, per, l in zip(phase, period, labels): ax.plot(per, a, label=l)

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"{vessel}/Load_RAO/f_ph_{i + 1}_head_{h}.png")
        plt.close(fig)



debug = True