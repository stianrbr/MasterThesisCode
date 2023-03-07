from misc import WADAM_res
import matplotlib.pyplot as plt

file = "S:\\Master\\Second order body convergence\\INO\\Mesh_4_m\\HydroDActivity1\\Analyses\\WadamAnalysis1\\WADAM1.LIS"


results1 = WADAM_res(file)

added_mass, period = results1.total_added_mass(2,2)



plt.plot(period, added_mass)
plt.show()

heave_amp, heave_ph, period = results1.motion_RAO(2,0)

plt.plot(period, heave_amp)
plt.show()

plt.plot(period, heave_ph)
plt.show()

debug = True