import numpy as np
import os
from Read_WADAM.misc import WADAM_res
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


def discretization_error(phi, phi_kappa, phi_kappa_squared, kappa):
    gamma = np.log10((phi_kappa-phi_kappa_squared)/(phi-phi_kappa))/np.log10(kappa)
    res = (phi-phi_kappa)/(np.float_power(np.full_like(phi, kappa), gamma)-1)
    return res





mainfolder = "S:\\Master\\Convergence study - body mesh\\"

vessel = "INO"

if vessel == "INO":
    title = "Peripheral tower design"
elif vessel == "Volturn":
    title = "Central tower design"
else:
    raise Exception("Check vessel")

kappa = 1.5

show = True





mesh_0_8 = mainfolder+vessel+"\\Mesh_0_8m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_0_8 = WADAM_res(filename=mesh_0_8)

mesh_1_2 = mainfolder+vessel+"\\Mesh_1_2m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_2 = WADAM_res(filename=mesh_1_2)

mesh_1_8 = mainfolder+vessel + "\\Mesh_1_8m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_1_8 = WADAM_res(filename=mesh_1_8)

mesh_2_7 = mainfolder +vessel+ "\\Mesh_2_7m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.lis"
res_2_7 = WADAM_res(filename=mesh_2_7)


results = [res_0_8, res_1_2, res_1_8, res_2_7]


try:
    os.mkdir("Relative_error")
except FileExistsError:
    pass

try:
    os.mkdir("Relative_error\\"+vessel)
except FileExistsError:
    pass


"""
Added mass relative error
"""

try:
    os.mkdir("Relative_error\\"+f"{vessel}\\Added_mass")
except FileExistsError:
    pass

a_m_pairs = [(0,0), (2,2), (4,4)]

colors = sns.color_palette("colorblind", n_colors=len(a_m_pairs))
fig = plt.figure()
ax = plt.axes()

markers = ["o", "v", "*"]
marker_flag = [False, False, False]
marker_size = 3

line_styles = ["solid", "dashed", "dotted"]

a_m_pair_label = [r"A$_{11}$", r"A$_{33}$", r"A$_{55}$"]

mesh_labels = ["0.8 m", "1.2 m", "1.8 m", "2.7 m"]

for pp, pair in enumerate(a_m_pairs):

    a_ij, per = zip(*[res.get_Total_added_mass(*pair) for res in results])
    rel_err_a_ij = []
    i=0
    while i+3 <= len(a_ij):
        rel_err_a_ij.append(discretization_error(*a_ij[i:i+3], kappa))
        i+=1

    for ii, arr in enumerate(rel_err_a_ij):
        ax.plot(per[ii], arr/a_ij[ii]*100, marker=markers[ii], markersize=marker_size, ls=line_styles[ii], c=colors[pp])


for lp, lab in enumerate(a_m_pair_label):
    ax.plot([], [], ls="solid", c=colors[lp], label=lab)

ax.plot([], [], ls="solid", c="white", label=" ")

for ii, arr in enumerate(rel_err_a_ij):
    ax.plot([], [], marker=markers[ii], ls=line_styles[ii], c="k", label=mesh_labels[ii])

plt.xlabel("Period [s]")
plt.ylabel("Relative error [\%]")
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.title(f"Convergence study - Added mass \n {title}")
plt.tight_layout()
plt.savefig("Relative_error\\"+f"{vessel}\\Added_mass\\a_rel_err.png")
if show:
    plt.show()
else:
    plt.close(fig)

"""
Load RAO relative error
"""

try:
    os.mkdir("Relative_error\\"+f"{vessel}\\Load_RAO")
except FileExistsError:
    pass

i_head_pairs = [(0,0), (0,2), (2,0), (4, 0)]

colors = sns.color_palette("colorblind", n_colors=len(i_head_pairs))
fig = plt.figure()
ax = plt.axes()

markers = ["o", "v", "*"]
marker_flag = [False, False, False]
marker_size = 3

line_styles = ["solid", "dashed", "dotted"]

i_head_pairs_label = [r"F$_{1}$ - 0$^{\circ}$", r"F$_{1}$ - 45$^{\circ}$", r"F$_{3}$ - 0$^{\circ}$", r"F$_{5}$ - 0$^{\circ}$"]

mesh_labels = ["0.8 m", "1.2 m", "1.8 m", "2.7 m"]

for pp, pair in enumerate(i_head_pairs):

    f_ij, ph_ij, per = zip(*[res.get_Load_RAO(*pair) for res in results])
    rel_err_f_ij = []
    i=0
    while i+3 <= len(f_ij):
        rel_err_f_ij.append(discretization_error(*f_ij[i:i+3], kappa))
        i+=1

    for ii, arr in enumerate(rel_err_f_ij):
        ax.plot(per[ii], arr/f_ij[ii]*100, marker=markers[ii], markersize=marker_size, ls=line_styles[ii], c=colors[pp])


for lp, lab in enumerate(i_head_pairs_label):
    ax.plot([], [], ls="solid", c=colors[lp], label=lab)

ax.plot([], [], ls="solid", c="white", label=" ")

for ii, arr in enumerate(rel_err_a_ij):
    ax.plot([], [], marker=markers[ii], ls=line_styles[ii], c="k", label=mesh_labels[ii])

plt.xlabel("Period [s]")
plt.ylabel("Relative error [\%]")
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.title(f"Convergence study - Wave load excitation \n {title}")
plt.tight_layout()
plt.savefig("Relative_error\\"+f"{vessel}\\Load_RAO\\f_rel_err.png")
if show:
    plt.show()
else:
    plt.close(fig)


"""
Motion RAO relative error
"""

try:
    os.mkdir("Relative_error\\"+f"{vessel}\\Motion_RAO")
except FileExistsError:
    pass

i_head_pairs = [(0,0), (0,2), (2,0), (4, 0)]

colors = sns.color_palette("colorblind", n_colors=len(i_head_pairs))
fig = plt.figure()
ax = plt.axes()

markers = ["o", "v", "*"]
marker_flag = [False, False, False]
marker_size = 3

line_styles = ["solid", "dashed", "dotted"]

i_head_pairs_label = [r"X$_{1}$ - 0$^{\circ}$", r"X$_{1}$ - 45$^{\circ}$", r"X$_{3}$ - 0$^{\circ}$", r"X$_{5}$ - 0$^{\circ}$"]

mesh_labels = ["0.8 m", "1.2 m", "1.8 m", "2.7 m"]

for pp, pair in enumerate(i_head_pairs):

    x_ij, ph_ij, per = zip(*[res.get_Motion_RAO(*pair) for res in results])
    rel_err_x_ij = []
    i=0
    while i+3 <= len(x_ij):
        rel_err_x_ij.append(discretization_error(*x_ij[i:i+3], kappa))
        i+=1

    for ii, arr in enumerate(rel_err_x_ij):
        ax.plot(per[ii], arr/x_ij[ii]*100, marker=markers[ii], markersize=marker_size, ls=line_styles[ii], c=colors[pp])


for lp, lab in enumerate(i_head_pairs_label):
    ax.plot([], [], ls="solid", c=colors[lp], label=lab)

ax.plot([], [], ls="solid", c="white", label=" ")

for ii, arr in enumerate(rel_err_a_ij):
    ax.plot([], [], marker=markers[ii], ls=line_styles[ii], c="k", label=mesh_labels[ii])

plt.xlabel("Period [s]")
plt.ylabel("Relative error [\%]")
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.title(f"Convergence study - Wave load excitation \n {title}")
plt.tight_layout()
plt.savefig("Relative_error\\"+f"{vessel}\\Motion_RAO\\x_rel_err.png")
if show:
    plt.show()
else:
    plt.close(fig)
