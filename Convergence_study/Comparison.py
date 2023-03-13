from Read_WADAM.misc import WADAM_res
from Read_WADAM.misc import Addedmassunits, Dampingunits, Excitation_units, Mean_drift_units, motionRAO_units

import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


mainfolder = "S:\\Master\\Convergence study - body mesh\\"

vessel = "Volturn"

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

mesh_1_5 = mainfolder+ vessel + "\\Mesh_1_5m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_5 = WADAM_res(filename=mesh_1_5)

mesh_1_8 = mainfolder+ vessel + "\\Mesh_1_8m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_8 = WADAM_res(filename=mesh_1_8)

mesh_2_7 = mainfolder+ vessel + "\\Mesh_2_7m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_2_7 = WADAM_res(filename=mesh_2_7)

results = [res_0_8, res_1_2, res_1_5, res_1_8, res_2_7]
labels = ["0.8 m", "1.2 m", "1.5 m", "1.8 m", "2.7 m"]

res_2_7.get_Mean_drift_full(0,0)
try:
    os.mkdir("Comparison")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+vessel)
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+f"{vessel}\\Added_mass")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+f"{vessel}\\Damping")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+f"{vessel}\\Motion_RAO")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+f"{vessel}\\Load_RAO")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison\\"+f"{vessel}\\Mean_drift")
except FileExistsError:
    pass


colors = sns.color_palette("colorblind", len(results))


for i in range(6):
    for j in range(6):

        added_mass, period = zip(*[res.get_Total_added_mass(i,j) for res in results])
        fig = plt.figure()
        ax = plt.axes()
        for ci, (a,p ,l) in enumerate(zip(added_mass, period, labels)): ax.plot(p, a, c=colors[ci], label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nAdded mass - (i,j)=({i+1},{j+1})")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Added mass [" + Addedmassunits[i][j] +"]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Added_mass\\a_{i+1}{j+1}.png")
        plt.close(fig)

        damping, period = zip(*[res.get_Total_damping(i,j) for res in results])

        fig = plt.figure()
        ax = plt.axes()
        for ci, (b, p, l) in enumerate(zip(damping, period, labels)): ax.plot(p, b, c=colors[ci],  label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nDamping - (i,j)=({i + 1},{j + 1})")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Damping [" + Dampingunits[i][j] + "]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Damping\\b_{i + 1}{j + 1}.png")
        plt.close(fig)

headings = results[0].Environment_data.wave_headings

for j, h in enumerate(headings):
    for i in range(6):
        amp, phase, period = zip(*[res.get_Motion_RAO(i,j) for res in results])

        fig = plt.figure()
        ax = plt.axes()
        for ci, (a, per, l) in enumerate(zip(amp, period, labels)): ax.plot(per, a, c=colors[ci], label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nMotion RAO amplitude - i={i + 1}\tHeading: {h}"+r"$^{\circ}$ ")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Displacement ["+motionRAO_units[i]+"]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Motion_RAO\\x_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        for ci, (a, per, l) in enumerate(zip(phase, period, labels)): ax.plot(per, a, c=colors[ci], label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nMotion RAO phase - i={i + 1}\tHeading: {h}"+r"$^{\circ}$ ")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Phase [deg]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Motion_RAO\\x_ph_{i + 1}_head_{h}.png")
        plt.close(fig)

        amp, phase, period = zip(*[res.get_Load_RAO(i,j) for res in results])


        fig = plt.figure()
        ax = plt.axes()
        for ci, (a, per, l) in enumerate(zip(amp, period, labels)): ax.plot(per, a, c=colors[ci], label=l)
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nLoad RAO amplitude - i={i + 1}\tHeading: {h}"+r"$^{\circ}$ ")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Excitation ["+Excitation_units[i]+"]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Load_RAO\\f_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        for ci, (a, per, l) in enumerate(zip(phase, period, labels)): ax.plot(per, a, c=colors[ci], label=l)

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nLoad RAO phase - i={i + 1}\tHeading: {h}"+r"$^{\circ}$ ")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Phase [deg]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Load_RAO\\f_ph_{i + 1}_head_{h}.png")
        plt.close(fig)

        amp, period = zip(*[res.get_Mean_drift_full(i,j ) for res in results])

        fig = plt.figure()
        ax = plt.axes()
        for ci, (a, per, l) in enumerate(zip(amp, period, labels)): ax.plot(per, a, c=colors[ci], ls="solid", label=l)
        if i in [0, 1, 5]:
            amp_h, period_h = zip(*[res.get_Mean_drift_horizontal(i,j) for res in results])
            for ci, (ah, ph, l) in enumerate(zip(amp, period, labels)): ax.plot(ph, ah, c=colors[ci], ls="dotted")
            ax.plot([], [], ls="solid", c="white", label=" ")
            ax.plot([], [], ls="solid", c="k", label="PI+CS")
            ax.plot([], [], ls="dotted", c="k", label="Far-field")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"{title}\nMean drift - i={i + 1}\tHeading: {h}" + r"$^{\circ}$ ")
        plt.xlabel(r"Period [s]")
        plt.ylabel(r"Mean drift ["+Mean_drift_units[i]+"]")
        plt.tight_layout()
        plt.savefig(f"Comparison\\{vessel}\\Mean_drift\\md_{i + 1}_head_{h}.png")
        plt.close(fig)





