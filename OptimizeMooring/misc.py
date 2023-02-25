import numpy as np
import OrcFxAPI

class Applied_load:
    def __init__(self, x=0, y=0, z=0, Fx=0, Fy=0, Fz=0):
        self.X = x
        self.Y = y
        self.Z = z
        self.FX = Fx
        self.FY = Fy
        self.FZ = Fz

class initital_pos:
    def __init__(self, x=0, y=0, z=0):
        self.X = x
        self.Y = y
        self.Z = z

class Constrained_model:
    def __init__(self, lindx, model, path, vesselName, mooringName, loads, w, p):
        self.lindx = lindx
        self.Model = model
        self.Path = path
        self.VesselName = vesselName
        self.MooringName = mooringName
        self.Loads = loads
        self.W = w
        self.P = p
    def objective(self, param):
        """
        Objective function: Returning the cost of a single mooring line
        :param param: (list-type) Contains variables in optimization problem
        :return: Scaled cost of single mooring line
        """
        scaling = 10000  # To scale value for faster convergence
        r, *l_vals = param

        cost = np.sum(np.pi / 4 * (np.array(l_vals) * self.W * self.P))

        print(f"r: {round(r, 4)} m - L: {np.around(l_vals, 4)} m - Cost : {round(cost/scaling, 2)}")
        return cost/scaling

    def RatedTension(self, param):
        scaling = 1000   # To scale value for faster convergence
        r, *l_vals = param

        self.Model.Reset()
        update_Mooring(r, l_vals, self.lindx, self.Model)
        self.Apply_loads()
        try:
            self.Model.CalculateStatics()
        except:
            self.Model.SaveData(self.Path + "Failed_statics.dat")
            return -np.inf

        tension = self.Model[self.MooringName].StaticResult("End force", objectExtra=OrcFxAPI.oeEndA)
        print(f"Rated tension: {tension/scaling}")
        return tension/scaling

    def PreTension(self, param):
        scaling = 1000  # To scale value for faster convergence
        r, *l_vals = param

        self.Model.Reset()
        update_Mooring(r, l_vals, self.lindx, self.Model)
        vessel = self.Model[self.VesselName]
        vessel.IncludeAppliedLoads = "No"
        try:
            self.Model.CalculateStatics()
        except:
            self.Model.SaveData(self.Path + "Failed_statics.dat")
            return -np.inf

        tension = self.Model[self.MooringName].StaticResult("End force", objectExtra=OrcFxAPI.oeEndA)
        print(f"Pretension: {tension/scaling}")
        return tension/scaling

    def static_offset(self, param):
        scaling = 10  # To scale value for faster convergence
        self.Model.Reset()
        r, *l_vals = param

        update_Mooring(r, l_vals, self.lindx, self.Model)

        vessel = self.Model[self.VesselName]
        self.Apply_loads()

        try:
            self.Model.CalculateStatics()
        except:
            return -np.inf
        new_x = vessel.StaticResult("X")
        print(f"Offset: {np.abs(new_x)}")
        return np.abs(new_x)/scaling

    def Modes_period(self, param):
        self.Model.Reset()
        r, *l_vals = param
        update_Mooring(r, l_vals, self.lindx, self.Model)
        vessel = self.Model[self.VesselName]
        vessel.IncludeAppliedLoads = "No"
        try:
            self.Model.CalculateStatics()
        except:
            return -np.inf
        mode_period = OrcFxAPI.Modes(self.Model,
                                     OrcFxAPI.ModalAnalysisSpecification(calculateShapes=False, firstMode=1,
                                                                         lastMode=1,
                                                                         includeCoupledObjects=True)).period[0]
        print(f"Period: {mode_period}")
        return mode_period/100

    def Rated_touchdown(self, param):
        self.Model.Reset()
        r, *l_vals = param
        update_Mooring(r, l_vals, self.lindx, self.Model)

        vessel = self.Model[self.VesselName]
        self.Apply_loads()

        try:
            self.Model.CalculateStatics()
        except:
            return -np.inf
        line_length = self.Model[self.MooringName].StaticResult("Arc length", objectExtra=OrcFxAPI.oeEndB)
        td_arc=self.Model[self.MooringName].StaticResult("Arc length", objectExtra=OrcFxAPI.oeTouchdown)
        print(f"Rated touchdown: {td_arc/line_length}")
        return td_arc/line_length

    def Apply_loads(self):
        vessel = self.Model[self.VesselName]
        if len(self.Loads) > 0:
            vessel.IncludeAppliedLoads = "Yes"
            vessel.NumberOfLocalAppliedLoads = len(self.Loads)
            for li, load in enumerate(self.Loads):
                vessel.LocalAppliedLoadOriginX[li] = load.X
                vessel.LocalAppliedLoadOriginY[li] = load.Y
                vessel.LocalAppliedLoadOriginZ[li] = load.Z
                vessel.LocalAppliedForceX[li] = load.FX
                vessel.LocalAppliedForceY[li] = load.FY
                vessel.LocalAppliedForceZ[li] = load.FZ
        else:
            vessel.IncludeAppliedLoads = "No"


def update_Mooring(r, l_vals, l_index, model):
    """
    :param param: (Tuple) Single parameter for running scipy
                    r: Radial anchor position from origin
                    l: Unstreched mooring line length
    :return: SSE of mooring line stiffness
    """
    scaling = 100
    r = r*scaling
    l_vals = [l*scaling for l in l_vals]
    model["Mooring1"].EndBx, model["Mooring1"].EndBy = r * np.cos(np.deg2rad(0)), r * np.sin(np.deg2rad(0))
    model["Mooring2"].EndBx, model["Mooring2"].EndBy = r * np.cos(np.deg2rad(240)), r * np.sin(np.deg2rad(240))
    model["Mooring3"].EndBx, model["Mooring3"].EndBy = r * np.cos(np.deg2rad(120)), r * np.sin(np.deg2rad(120))

    for l, index in zip(l_vals, l_index):
        model["Mooring1"].length[index] = l
        model["Mooring2"].length[index] = l
        model["Mooring3"].length[index] = l
    print("Updated mooring")