from Read_WADAM.misc import WADAM_res
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


mainfolder = "S:\\Master\\Body mesh - sensitivity full solution\\INO\\"

mesh_0_5 = mainfolder+"Mesh_0_5_linear\\Mesh_0_5_linear\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_0_5 = WADAM_res(filename=mesh_0_5)

mesh_1 = mainfolder+"Mesh_1_0_linear\\Mesh_1_linear\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1 = WADAM_res(filename=mesh_1)

mesh_2 = mainfolder+"Mesh_2_full\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_2 = WADAM_res(filename=mesh_2)

mesh_4 =mainfolder+"Mesh_4_full\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_4 = WADAM_res(filename=mesh_4)

try:
    os.mkdir("INO")
except FileExistsError:
    pass

try:
    os.mkdir("INO\\Added_mass")
except FileExistsError:
    pass

try:
    os.mkdir("INO\\Damping")
except FileExistsError:
    pass

try:
    os.mkdir("INO\\Motion_RAO")
except FileExistsError:
    pass

try:
    os.mkdir("INO\\Load_RAO")
except FileExistsError:
    pass

for i in range(6):
    for j in range(6):
        a_0_5, p_0_5 = res_0_5.get_Total_added_mass(i, j)
        a_1, p_1 = res_1.get_Total_added_mass(i, j)
        a_2, p_2 = res_2.get_Total_added_mass(i, j)
        a_4, p_4 = res_4.get_Total_added_mass(i, j)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, a_0_5, label="0.5 m")
        ax.plot(p_1, a_1, label="1.0 m")
        ax.plot(p_2, a_2, label="2.0 m")
        ax.plot(p_4, a_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Added mass - (i,j)=({i+1},{j+1})")
        plt.tight_layout()
        plt.savefig(f"INO/Added_mass/a_{i+1}{j+1}.png")
        plt.close(fig)

        b_0_5, p_0_5 = res_0_5.get_Total_damping(i, j)
        b_1, p_1 = res_1.get_Total_damping(i, j)
        b_2, p_2 = res_2.get_Total_damping(i, j)
        b_4, p_4 = res_4.get_Total_damping(i, j)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, b_0_5, label="0.5 m")
        ax.plot(p_1, b_1, label="1.0 m")
        ax.plot(p_2, b_2, label="2.0 m")
        ax.plot(p_4, b_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Damping - (i,j)=({i + 1},{j + 1})")
        plt.tight_layout()
        plt.savefig(f"INO/Damping/b_{i + 1}{j + 1}.png")
        plt.close(fig)

headings = res_0_5.Environment_data.wave_headings

for j, h in enumerate(headings):
    for i in range(6):
        x_0_5, x_ph_0_5, p_0_5 = res_0_5.get_Motion_RAO(i, j)
        x_1, x_ph_1, p_1 = res_1.get_Motion_RAO(i, j)
        x_2, x_ph_2, p_2 = res_2.get_Motion_RAO(i, j)
        x_4, x_ph_4, p_4 = res_4.get_Motion_RAO(i, j)


        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, x_0_5, label="0.5 m")
        ax.plot(p_1, x_1, label="1.0 m")
        ax.plot(p_2, x_2, label="2.0 m")
        ax.plot(p_4, x_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"INO/Motion_RAO/x_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, x_ph_0_5, label="0.5 m")
        ax.plot(p_1, x_ph_1, label="1.0 m")
        ax.plot(p_2, x_ph_2, label="2.0 m")
        ax.plot(p_4, x_ph_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"INO/Motion_RAO/x_ph_{i + 1}_head_{h}.png")
        plt.close(fig)

        f_0_5, f_ph_0_5, p_0_5 = res_0_5.get_Load_RAO(i, j)
        f_1, f_ph_1, p_1 = res_1.get_Load_RAO(i, j)
        f_2, f_ph_2, p_2 = res_2.get_Load_RAO(i, j)
        f_4, f_ph_4, p_4 = res_4.get_Load_RAO(i, j)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, f_0_5, label="0.5 m")
        ax.plot(p_1, f_1, label="1.0 m")
        ax.plot(p_2, f_2, label="2.0 m")
        ax.plot(p_4, f_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"INO/Load_RAO/f_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_0_5, x_ph_0_5, label="0.5 m")
        ax.plot(p_1, f_ph_1, label="1.0 m")
        ax.plot(p_2, f_ph_2, label="2.0 m")
        ax.plot(p_4, f_ph_4, label="4.0 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"INO/Load_RAO/f_ph_{i + 1}_head_{h}.png")
        plt.close(fig)



debug = True