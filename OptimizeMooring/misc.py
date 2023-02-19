import numpy as np

def update_Mooring(param, args):
    """
    :param param: (Tuple) Single parameter for running scipy
                    r: Radial anchor position from origin
                    l: Unstreched mooring line length
    :return: SSE of mooring line stiffness
    """
    l = param

    new_model, new_path, ref, r = args

    new_vessel = new_model["VolturnUS-S"]


    new_model["Mooring1"].length[0] = l
    new_model["Mooring2"].length[0] = l
    new_model["Mooring3"].length[0] = l

    new_vessel.IncludeAppliedLoads = "Yes"

    new_vessel.NumberOfLocalAppliedLoads = 1

    new_vessel.LocalAppliedLoadOriginZ[0] = 150

    new_vessel.LocalAppliedForceX[0] = -2500  # [kN] rated thrust

    try:
        new_model.CalculateStatics()
    except:
        new_model.SaveData(new_path+"Failed_statics.dat")

    new_x = new_vessel.StaticResult("X")

    SSE = np.abs(new_x-ref)
    print(f"Radius: {r} - L: {l} - Error: {SSE}")
    return SSE