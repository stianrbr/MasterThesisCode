import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import OrcFxAPI
from misc import update_Mooring, Applied_load, Constrained_model

"""
Loading reference results
"""


ref_path = "S:\\Master\\Mooring_system\\Original_model\\"
ref_file = "K03 15MW semi-sub FOWT.dat"
ref_model = OrcFxAPI.Model(ref_path+ref_file)
ref_model.CalculateStatics()
ref_vessel = ref_model["VolturnUS-S"]
ref_end_force = ref_model["Mooring1"].StaticResult("End force", objectExtra=OrcFxAPI.oeEndA)
ref_modes = OrcFxAPI.Modes(ref_model, OrcFxAPI.ModalAnalysisSpecification(calculateShapes=False, firstMode=1, lastMode=1, includeCoupledObjects=True))

"""
Loading new model
"""


new_path = "S:\\Master\\Mooring_system\\INO moor\\"
new_file = "K03 15MW semi-sub FOWT.dat"
new_model = OrcFxAPI.Model(new_path+new_file)
new_vessel = new_model["VolturnUS-S"]
new_model.environment.WaterDepth = 300.0


l_index = (2,3)

optimization = Constrained_model(lindx=l_index,
                                 model=new_model,
                                 path=new_path,
                                 vesselName="VolturnUS-S",
                                 mooringName="Mooring1",
                                 loads=[Applied_load(z=150, Fx=-1790)],
                                 w=np.array([0.2825, 3299.2]),
                                 p=np.array([7.00, 2.50]))

"""
Initial guess, bounds and constraints are scaled to obtain values of the order of E+01-E+02,
apparently this improves the convergence of the SLSQP-method
"""

initital_guess = [1020/100, 345.0/100, 579.0/100]

bounds = scipy.optimize.Bounds([900.0/100, 100.0/100, 400.0/100], [1200.0/100, 400.0/100, 1000.0/100], keep_feasible=True)

constraints = [scipy.optimize.NonlinearConstraint(optimization.static_offset, 0.0/10, 18.0/10, keep_feasible=True),
               scipy.optimize.NonlinearConstraint(optimization.Modes_period, 100.0/100, 150.0/100, keep_feasible=True),
               scipy.optimize.NonlinearConstraint(optimization.PreTension, 1000.0/1000, 3500.0/1000, keep_feasible=True),
               scipy.optimize.NonlinearConstraint(optimization.RatedTension, 0.0/1000, 5000.0/1000, keep_feasible=True),
               scipy.optimize.NonlinearConstraint(optimization.Rated_touchdown, 0.0/100, 80.0/100, keep_feasible=True)]

#result = scipy.optimize.minimize(optimization.objective, initital_guess, method='trust-constr', bounds=bounds, constraints=constraints, options= {"maxiter":2000,
#                                                                                                                                                  "xtol":0.0001,
#                                                                                                                                                  "barrier_tol":0.0001,
#                                                                                                                                                  'verbose': 3,
#                                                                                                                                                  'initial_constr_penalty':100})

result = scipy.optimize.minimize(optimization.objective, initital_guess, method='SLSQP', bounds=bounds, constraints=constraints, options= {"maxiter":2000,
                                                                                                                                                  'verbose': 3})

if result.success:
    print(result.x)
    print(optimization.objective(result.x))
    update_Mooring(result.x[0], (result.x[1], result.x[2]), optimization.lindx, optimization.Model)
    optimization.Model.SaveData(new_path+"Optimized_newlines.dat")
else:
    print(result.message)
debug = True