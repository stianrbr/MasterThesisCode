import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize
from misc import update_Mooring
import OrcFxAPI



"""
Loading reference results
"""

ref_path = "S:\\Master\\Mooring_system\\Original_model\\"

ref_file = "K03 15MW semi-sub FOWT.dat"

ref_model = OrcFxAPI.Model(ref_path+ref_file)


ref_vessel = ref_model["VolturnUS-S"]

ref_vessel.IncludeAppliedLoads = "Yes"

ref_vessel.NumberOfLocalAppliedLoads = 1

ref_vessel.LocalAppliedLoadOriginZ[0] = 150

ref_vessel.LocalAppliedForceX[0] = -2500  # [kN] rated thrust

ref_model.CalculateStatics()

ref_x = ref_vessel.StaticResult("X")


"""
Loading new model
"""

new_path = "S:\\Master\\Mooring_system\\New_model\\"

new_file = "K03 15MW semi-sub FOWT.dat"

new_model = OrcFxAPI.Model(new_path+new_file)



new_model.environment.WaterDepth = 300.0

res = []
res_tension = []
df = pd.DataFrame(columns=["r", "l", "Pretension"])
r_arr = np.arange(800, 1200, 10)

for r in r_arr:
    new_model["Mooring1"].EndBx, new_model["Mooring1"].EndBy = r * np.cos(np.deg2rad(0)), r * np.sin(np.deg2rad(0))
    new_model["Mooring2"].EndBx, new_model["Mooring2"].EndBy = r * np.cos(np.deg2rad(240)), r * np.sin(np.deg2rad(240))
    new_model["Mooring3"].EndBx, new_model["Mooring3"].EndBy = r * np.cos(np.deg2rad(120)), r * np.sin(np.deg2rad(120))

    initial_guess = [r*1.05]

    args = ([new_model, new_path, ref_x, r])

    try:
        result = scipy.optimize.minimize(update_Mooring, initial_guess, args= args, method="Nelder-Mead", bounds=[(int(0.9*r),1200)], options= {"maxiter":1000, "xatol":0.01})
        if result.success:
            print(result.x)
            new_model.SaveData(new_path + "Optimized.dat")
            new_vessel = new_model["VolturnUS-S"]
            new_vessel.IncludeAppliedLoads = "No"
            new_model.CalculateStatics()
            pretension = new_model["Mooring1"].StaticResult("End force", objectExtra=OrcFxAPI.oeEndA)

            temp = {"r": [r], "l": [result.x[0]], "Pretension": [pretension]}
            df = pd.concat([df, pd.DataFrame(temp)])
        else:
                print("Failed")
    except:
        print("Failed")
        pass


df.to_hdf("Optimize_mooring.h5", key="df")