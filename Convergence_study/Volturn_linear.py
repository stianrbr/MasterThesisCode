from Read_WADAM.misc import WADAM_res
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


mainfolder = "S:\\Master\\Convergence study - body mesh\\Volturn\\"

start_period = 4
end_period = 50


mesh_1_5 = mainfolder+"Mesh_1_5m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_5 = WADAM_res(filename=mesh_1_5)

start_slice = np.argwhere(res_1_5.Environment_data.wave_periods>=start_period)[0][0]
end_slice = np.argwhere(res_1_5.Environment_data.wave_periods<=end_period)[-1][0]

start_slice = 0
end_slice = -1

mesh_2_25 = mainfolder + "Mesh_2_25m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_2_25 = WADAM_res(filename=mesh_2_25)

mesh_3_375 = mainfolder + "Mesh_3_375m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_3_375 = WADAM_res(filename=mesh_3_375)

try:
    os.mkdir("Comparison/Volturn")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/Volturn\\Added_mass")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/Volturn\\Damping")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/Volturn\\Motion_RAO")
except FileExistsError:
    pass

try:
    os.mkdir("Comparison/Volturn\\Load_RAO")
except FileExistsError:
    pass

for i in range(6):
    for j in range(6):
        a_1_5, p_1_5 = res_1_5.get_Total_added_mass(i, j)
        a_2_25, p_2_25 = res_2_25.get_Total_added_mass(i, j)
        a_3_375, p_3_375 = res_3_375.get_Total_added_mass(i, j)


        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], a_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], a_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], a_3_375[start_slice:end_slice], label="3.375 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Added mass - (i,j)=({i+1},{j+1})")
        plt.tight_layout()
        plt.savefig(f"Volturn/Added_mass/a_{i+1}{j+1}.png")
        plt.close(fig)

        b_1_5, p_1_5 = res_1_5.get_Total_damping(i, j)
        b_2_25, p_2_25 = res_2_25.get_Total_damping(i, j)
        b_3_375, p_3_375 = res_3_375.get_Total_damping(i, j)


        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], b_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], b_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], b_3_375[start_slice:end_slice], label="3.375 m")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Damping - (i,j)=({i + 1},{j + 1})")
        plt.tight_layout()
        plt.savefig(f"Volturn/Damping/b_{i + 1}{j + 1}.png")
        plt.close(fig)

headings = res_1_5.Environment_data.wave_headings

for j, h in enumerate(headings):
    for i in range(6):
        x_1_5, x_ph_1_5, p_1_5 = res_1_5.get_Motion_RAO(i, j)
        x_2_25, x_ph_2_25, p_2_25 = res_2_25.get_Motion_RAO(i, j)
        x_3_375, x_ph_3_375, p_3_375 = res_3_375.get_Motion_RAO(i, j)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], x_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], x_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], x_3_375[start_slice:end_slice], label="3.375 m")

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"Volturn/Motion_RAO/x_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], x_ph_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], x_ph_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], x_ph_3_375[start_slice:end_slice], label="3.375 m")

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Motion RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"Volturn/Motion_RAO/x_ph_{i + 1}_head_{h}.png")
        plt.close(fig)

        f_1_5, f_ph_1_5, p_1_5 = res_1_5.get_Load_RAO(i, j)
        f_2_25, f_ph_2_25, p_2_25 = res_2_25.get_Load_RAO(i, j)
        f_3_375, f_ph_3_375, p_3_375 = res_3_375.get_Load_RAO(i, j)


        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], f_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], f_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], f_3_375[start_slice:end_slice], label="3.375 m")

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO amplitude - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"Volturn/Load_RAO/f_{i + 1}_head_{h}.png")
        plt.close(fig)

        fig = plt.figure()
        ax = plt.axes()
        ax.plot(p_1_5[start_slice:end_slice], x_ph_1_5[start_slice:end_slice], label="1.5 m")
        ax.plot(p_2_25[start_slice:end_slice], f_ph_2_25[start_slice:end_slice], label="2.25 m")
        ax.plot(p_3_375[start_slice:end_slice], f_ph_3_375[start_slice:end_slice], label="3.375 m")

        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.title(f"Load RAO phase - i={i + 1}\nHeading: {h}"+r"$^{\circ}$ ")
        plt.tight_layout()
        plt.savefig(f"Volturn/Load_RAO/f_ph_{i + 1}_head_{h}.png")
        plt.close(fig)



debug = True