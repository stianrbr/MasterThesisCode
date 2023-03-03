import re
import numpy as np
import warnings
from io import StringIO
regexp = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"

class Radiation_res:
    def __init__(self,matrix=None, period=None):
        self.Matrix = matrix
        self.period = period

class Dimensions:
    def __init__(self, rho=None, g=None, vol=None, l=None, wa=None):
        self.rho = rho
        self.g = g
        self.vol = vol
        self.l  = l
        self.wa = wa

class Point:
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z

class Mass_matrix:
    def __init__(self, M, XG, YG, ZG, XRAD, YRAD, ZRAD, XYRAD, XZRAD, YZRAD):
        self.M = M
        self.XG, self.YG, self.ZG = XG, YG, ZG
        self.COG = Point(XRAD, YRAD, ZRAD)
        self.XYRAD, self.XZRAD, self.YZRAD = XYRAD, XZRAD, YZRAD


class Hydro_static:
    def __init__(self, vol_disp, t_moor_vertical, wp_area, xb, yb, zb):
        self.vol_disp = vol_disp
        self.t_moor_vertical = t_moor_vertical
        self.wp_area = wp_area
        self.COB = Point(xb, yb, zb)

class Environmental_data:
    def __init__(self, water_depth, wave_periods, wave_length):
        self.water_depth = water_depth
        self.wave_periods = wave_periods
        self.wave_lengths = wave_length

class Tank:
    def __init__(self, tank_numb=None, mass=None, cog=None, added_mass=None, damping=None, restoring=None):
        self.tank_numb = tank_numb
        self.mass = mass
        self.cog = cog
        self.added_mass = added_mass
        self.damping = damping
        self.restoring = restoring


def extract_number_from_lines(file, line_indx, numb_indx):
    temp = [float(re.findall(regexp, file[i])[j]) for i, j in zip(line_indx, numb_indx)]
    return temp
class WADAM_res:
    def __init__(self, filename):
        self.filename = filename
        self.compartment_flag = False
        self.mooring_flag = False
        # Added mass
        self.Total_added_mass = []
        self.Hull_added_mass = None
        self.Zero_freq_added_mass = None
        self.Inf_freq_added_mass = None
        # Damping
        self.Total_damping = []
        self.Hull_damping = None

        # Restoring
        self.Hydrostat_restoring = None
        self.Mooring_restoring = None
        self.Total_restoring = None
        # Load RAO
        self.loadRAO = None
        # Response RAO
        self.motionRAO = None
        # Tanks
        self.tanks = None

    def read_file(self):
        with open(self.filename) as f:
            self.datafile = f.readlines()

            for i, line in enumerate(self.datafile):
                if not line.isspace():
                    if "COMPARTMENT" in line:
                        self.compartment_flag = True
                        self.tanks = []
                        self.Hull_added_mass = []
                        self.Hull_damping = []
                    if "MOORING" in line:
                        self.mooring_flag = True
                    if "MASS PROPERTIES AND STRUCTURAL DATA:" in line:
                        try:
                            section = self.datafile[i:i+20]
                            m = float(re.findall(regexp,[s for s in section if "MASS OF THE STRUCTURE                      M" in s][0])[0])
                            xg = float(re.findall(regexp,[s for s in section if "CENTRE OF GRAVITY                          XG" in s][0])[0])
                            yg = float(re.findall(regexp, [s for s in section if "                                           YG" in s][0])[0])
                            zg = float(re.findall(regexp, [s for s in section if "                                           ZG" in s][0])[0])
                            xrad = float(re.findall(regexp, [s for s in section if "ROLL  RADIUS OF GYRATION                   XRAD" in s][0])[0])
                            yrad = float(re.findall(regexp, [s for s in section if "PITCH RADIUS OF GYRATION                   YRAD" in s][0])[0])
                            zrad = float(re.findall(regexp, [s for s in section if "YAW   RADIUS OF GYRATION                   ZRAD" in s][0])[0])
                            xyrad = float(re.findall(regexp, [s for s in section if "ROLL-PITCH CENTRIFUGAL MOMENT              XYRAD" in s][0])[0])
                            xzrad = float(re.findall(regexp, [s for s in section if "ROLL-YAW   CENTRIFUGAL MOMENT              XZRAD" in s][0])[0])
                            yzrad = float(re.findall(regexp, [s for s in section if "PITCH-YAW  CENTRIFUGAL MOMENT              YZRAD" in s][0])[0])
                            self.Massmatrix = Mass_matrix(m, xg, yg, zg, xrad, yrad, zrad, xyrad, xzrad, yzrad)
                        except:
                            raise Exception("Failed reading mass properties")
                    elif "HYDROSTATIC DATA:" in line:
                        try:
                            section = self.datafile[i:i + 20]
                            vol_disp = float(re.findall(regexp,[s for s in section if "DISPLACED VOLUME                           VOL" in s][0])[0])
                            try:
                                t_moor_vertical = float(re.findall(regexp,[s for s in section if "VERTICAL PRETENSION OF MOORING ELEMENTS" in s][0])[0])
                            except:
                                t_moor_vertical = None
                                warnings.warn("No vertical mooring tension found.")
                            wp_area = float(re.findall(regexp,[s for s in section if "WATER PLANE AREA                           WPLA" in s][0])[0])
                            xb = float(re.findall(regexp,[s for s in section if "CENTRE OF BUOYANCY                         XCB" in s][0])[0])
                            yb = float(re.findall(regexp,[s for s in section if "                                           YCB" in s][0])[0])
                            zb = float(re.findall(regexp,[s for s in section if "                                           ZCB" in s][0])[0])
                            self.Hydro_static = Hydro_static(vol_disp, wp_area, t_moor_vertical, xb, yb, zb)
                        except:
                            raise Exception("Failed reading hydrostatic properties")
                    elif "2.8 ENVIRONMENTAL DATA:" in line:
                        try:
                            section = self.datafile[i:i+10]
                            water_depth = float(re.findall(regexp,[s for s in section if "WATER DEPTH                                          =" in s][0])[0])
                            n_wave_per = int(re.findall(regexp,[s for s in section if "NUMBER OF WAVE LENGTHS                               =" in s][0])[0])
                            n_wave_head = int(re.findall(regexp,[s for s in section if "NUMBER OF HEADING ANGLES                             =" in s][0])[0])

                        except:
                            raise Exception("Failed reading environmental data")
                    elif "WAVE DESCRIPTION:" in line:
                        try:
                            section = self.datafile[i+5:i+n_wave_per+5]
                            wave_per = [float(re.findall(regexp, w)[3]) for w in section]
                        except:
                            raise Exception("Failed reading wave periods")
                    elif "HEADING ANGLES (ANGLE BETWEEN POS. X-AXIS AND DIRECTION" in line:
                        try:
                            section = self.datafile[i+5:i+n_wave_head+5]
                            wave_head = [float(re.findall(regexp, w)[1]) for w in section]
                            self.Environment_data = Environmental_data(water_depth, np.array(wave_per), np.array(wave_head))
                        except:
                            raise Exception("Failed reading wave headings")

                    """
                    Extracting restoring data
                    """

                    if "HYDROSTATIC RESTORING COEFFICIENT MATRIX" in line:
                        c_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(float)
                        self.Hydrostat_restoring = c_matrix
                    if "TOTAL RESTORING COEFFICIENT MATRIX" in line:
                        c_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(float)
                        self.Total_restoring = c_matrix
                    if self.compartment_flag:
                        if "STATIC RESTORING COEFFICIENT MATRIX FOR FLUID IN TANK" in line:
                            tank_numb = int(re.findall(regexp, line)[0])
                            c_t_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(float)
                            self.tanks.append(Tank(tank_numb=tank_numb, restoring=c_t_matrix))
                    if self.mooring_flag:
                        if "STATIC RESTORING COEFFICIENT MATRIX FOR TLP AND MOORING ELEMENTS" in line:
                            self.Mooring_restoring = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(float)

                    """
                    Extracting added mass data
                    """

                    if "TOTAL ADDED MASS MATRIX " in line:
                        a_matrix = np.array([re.findall(regexp, self.datafile[i+4+a_i])[1:] for a_i in range(0,6)]).astype(float)
                        wave_len_number = int(re.findall(regexp, self.datafile[i-9])[1])
                        self.Total_added_mass.append(Radiation_res(matrix=a_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))

                    if "POTENTIAL ADDED MASS MATRIX FOR INFINITE WAVE PERIOD" in line:
                        self.Inf_freq_added_mass = np.array([re.findall(regexp, self.datafile[i+4+a_i])[1:] for a_i in range(0,6)]).astype(float)

                    if "POTENTIAL ADDED MASS MATRIX FOR ZERO WAVE PERIOD" in line:
                        self.Zero_freq_added_mass = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(float)

                    if self.compartment_flag:
                        if "ADDED MASS MATRIX OF OUTER HULL" in line:
                            a_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(float)
                            wave_len_number = int(re.findall(regexp, self.datafile[i - 9])[1])
                            self.Hull_added_mass.append(Radiation_res(matrix=a_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))
                        if "ADDED MASS MATRIX OF TANK " in line:
                            tank_numb = int(re.findall(regexp, line)[0])
                            a_t_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(float)
                            wave_len_number = int(re.findall(regexp, self.datafile[i - 9])[1])
                            assert(self.tanks[tank_numb-1].tank_numb == tank_numb)
                            self.tanks[tank_numb-1].added_mass = Radiation_res(matrix=a_t_matrix, period=self.Environment_data.wave_periods[wave_len_number-1])

                    """
                    Extracting damping data
                    """











